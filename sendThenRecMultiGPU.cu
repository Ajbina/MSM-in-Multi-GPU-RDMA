#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <arpa/inet.h>
#include <infiniband/verbs.h>
#include <unistd.h>
#include <string.h>
#include <cuda_runtime.h>

#define TCP_PORT 18515
#define RECEIVER1_IP "128.105.144.57"
#define RECEIVER2_IP "128.105.144.47"
#define RECEIVER3_IP "128.105.144.43"
#define RECEIVER4_IP "128.105.144.53"

bool pointInVector(const std::vector<ECPoint>& points, int x, int y) {
    for (const auto& point : points) {
        if (point.x == x && point.y == y) {
            return true;
        }
    }
    return false;
}

int modCPU(int a, int m) {
    int r = a % m;
    return (r < 0) ? r + m : r;
} 

int modInverseCPU(int a, int p) {
    a = modCPU(a, p);
    for (int x = 1; x < p; x++) {
        if (modCPU(a * x, p) == 1) {
            return x;
        }
    }
    return 0;
}

ECPoint ecDoubleTrueCPU(const ECPoint &P, int a, int p) {
    ECPoint R;
    if (P.y == 0){
        return {-1, -1, 1};
    }
    int infP = P.infinityFlag;
    if (infP == 1){
        return {-1, -1, 1};
    }

    int numerator = modCPU(3 * P.x * P.x + a, p);
    int denominator = modCPU(2 * P.y, p);
    int invDenom = modInverseCPU(denominator, p);
    int lambda = modCPU(numerator * invDenom, p);

    R.x = modCPU(lambda * lambda - 2 * P.x, p);
    R.y = modCPU(lambda * (P.x - R.x) - P.y, p);
    R.infinityFlag = 0;
    return R;
}

ECPoint ecAddTrueCPU(const ECPoint &P, const ECPoint &Q, int a, int p) {
    ECPoint R;
    int infP = P.infinityFlag;
    int infQ = Q.infinityFlag;
    if (infP == 1 && infQ == 1){
        return {-1, -1, 1};
    }
    if (infP == 1){
        return Q;
    }
    if (infQ == 1){
        return P;
    }
    if (P.x == Q.x && P.y != Q.y){
        return {-1, -1, 1};
    }
    if(P.x == Q.x && P.y == Q.y){
        R = ecDoubleTrueCPU(P, a, p);
        return R;
    }

    int numerator = modCPU(Q.y - P.y, p);
    int denominator = modCPU(Q.x - P.x, p);
    int invDenom = modInverseCPU(denominator, p);
    int lambda = modCPU(numerator * invDenom, p);

    R.x = modCPU(lambda * lambda - P.x - Q.x, p);
    R.y = modCPU(lambda * (P.x - R.x) - P.y, p);
    R.infinityFlag = 0;
    return R;
}

struct rdma_info {
    uint64_t addr;
    uint32_t rkey;
    uint32_t qp_num;
    uint16_t lid;
    uint8_t gid[16];
};

struct ECPoint {
    int x;
    int y;
    int infinityFlag;
};

struct RDMAConnection {
    ibv_qp* qp;
    ibv_mr* mr;
    ECPoint* buffer;
    rdma_info remote;
    int sockfd;
};

void modify_qp_to_rts(ibv_qp *qp, uint32_t remote_qpn, uint16_t remote_lid, uint8_t *remote_gid) {
    ibv_qp_attr attr = {};
   
    // Change QP state to INIT
    attr.qp_state = IBV_QPS_INIT;
    attr.pkey_index = 0;
    attr.port_num = 1;
    attr.qp_access_flags = IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE;
    ibv_modify_qp(qp, &attr, IBV_QP_STATE | IBV_QP_PKEY_INDEX | IBV_QP_PORT | IBV_QP_ACCESS_FLAGS);
   
    // Change QP state to RTR
    attr.qp_state = IBV_QPS_RTR;
    attr.path_mtu = IBV_MTU_1024;
    attr.dest_qp_num = remote_qpn;
    attr.rq_psn = 0;
    attr.max_dest_rd_atomic = 1;
    attr.min_rnr_timer = 12;
    attr.ah_attr.is_global = 1;
    attr.ah_attr.dlid = remote_lid;
    attr.ah_attr.sl = 0;
    attr.ah_attr.src_path_bits = 0;
    attr.ah_attr.port_num = 1;
    memcpy(attr.ah_attr.grh.dgid.raw, remote_gid, 16);
    attr.ah_attr.grh.flow_label = 0;
    attr.ah_attr.grh.hop_limit = 1;
    attr.ah_attr.grh.sgid_index = 0;
    attr.ah_attr.grh.traffic_class = 0;
   
    ibv_modify_qp(qp, &attr, IBV_QP_STATE | IBV_QP_AV | IBV_QP_PATH_MTU | IBV_QP_DEST_QPN |
                  IBV_QP_RQ_PSN | IBV_QP_MAX_DEST_RD_ATOMIC | IBV_QP_MIN_RNR_TIMER);
   
    // Change QP state to RTS
    attr.qp_state = IBV_QPS_RTS;
    attr.timeout = 14;
    attr.retry_cnt = 7;
    attr.rnr_retry = 7;
    attr.sq_psn = 0;
    attr.max_rd_atomic = 1;
   
    ibv_modify_qp(qp, &attr, IBV_QP_STATE | IBV_QP_TIMEOUT | IBV_QP_RETRY_CNT |
                  IBV_QP_RNR_RETRY | IBV_QP_SQ_PSN | IBV_QP_MAX_QP_RD_ATOMIC);
}

RDMAConnection setup_connection(const std::string& receiver_ip, ECPoint* data_chunk, size_t chunk_bytes, ibv_context* ctx, ibv_pd* pd, ibv_cq* cq) {
    RDMAConnection conn;
    conn.buffer = data_chunk;
    conn.mr = ibv_reg_mr(pd, data_chunk, chunk_bytes, IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_ZERO_BASED);

    ibv_qp_init_attr qp_init_attr = {};
    qp_init_attr.send_cq = cq;
    qp_init_attr.recv_cq = cq;
    qp_init_attr.qp_type = IBV_QPT_RC;
    qp_init_attr.cap.max_send_wr = 1;
    qp_init_attr.cap.max_recv_wr = 1;
    qp_init_attr.cap.max_send_sge = 1;
    qp_init_attr.cap.max_recv_sge = 1;

    conn.qp = ibv_create_qp(pd, &qp_init_attr);

    ibv_port_attr port_attr;
    ibv_query_port(ctx, 1, &port_attr);
    uint16_t local_lid = port_attr.lid;
    union ibv_gid gid;
    ibv_query_gid(ctx, 1, 0, &gid);

    // TCP exchange
    int sockfd = socket(AF_INET, SOCK_STREAM, 0);
    sockaddr_in servaddr{};
    servaddr.sin_family = AF_INET;
    servaddr.sin_port = htons(TCP_PORT);
    inet_pton(AF_INET, receiver_ip.c_str(), &servaddr.sin_addr);
    connect(sockfd, (sockaddr*)&servaddr, sizeof(servaddr));

    rdma_info local_info = {
        .addr = (uint64_t)(uintptr_t)data_chunk,
        .rkey = conn.mr->rkey,
        .qp_num = conn.qp->qp_num,
        .lid = local_lid
    };
    memcpy(local_info.gid, gid.raw, 16);
    send(sockfd, &local_info, sizeof(local_info), 0);
    recv(sockfd, &conn.remote, sizeof(conn.remote), 0);
    //close(sockfd);

    // Modify QP state
    modify_qp_to_rts(conn.qp, conn.remote.qp_num, conn.remote.lid, conn.remote.gid);

    conn.sockfd = sockfd;
    return conn;
}

int main() {
    std::ifstream file("data/points.csv");
    std::string line;
    std::getline(file, line);
    int numPointsCSV = std::stoi(line);

    std::vector<ECPoint> h_points;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<int> coords;
        while (std::getline(ss, token, ',')) coords.push_back(std::stoi(token));
        h_points.push_back({coords[0], coords[1], coords[2]});
    }
    file.close();

    int total = h_points.size();
    int quarter = total / 4;

    ibv_device** dev_list = ibv_get_device_list(nullptr);
    ibv_context* ctx = ibv_open_device(dev_list[0]);
    ibv_pd* pd = ibv_alloc_pd(ctx);
    ibv_cq* cq = ibv_create_cq(ctx, 10, nullptr, nullptr, 0);

    ECPoint* chunk1 = static_cast<ECPoint*>(malloc(quarter * sizeof(ECPoint)));
    ECPoint* chunk2 = static_cast<ECPoint*>(malloc(quarter * sizeof(ECPoint)));
    ECPoint* chunk3 = static_cast<ECPoint*>(malloc(quarter * sizeof(ECPoint)));
    ECPoint* chunk4 = static_cast<ECPoint*>(malloc((total - 3*quarter) * sizeof(ECPoint)));
    
    memcpy(chunk1, h_points.data(), quarter * sizeof(ECPoint));
    memcpy(chunk2, h_points.data() + quarter, quarter * sizeof(ECPoint));
    memcpy(chunk3, h_points.data() + 2*quarter, quarter * sizeof(ECPoint));
    memcpy(chunk4, h_points.data() + 3*quarter, (total - 3*quarter) * sizeof(ECPoint));

    RDMAConnection conn1 = setup_connection(RECEIVER1_IP, chunk1, quarter * sizeof(ECPoint), ctx, pd, cq);
    RDMAConnection conn2 = setup_connection(RECEIVER2_IP, chunk2, quarter * sizeof(ECPoint), ctx, pd, cq);
    RDMAConnection conn3 = setup_connection(RECEIVER3_IP, chunk3, quarter * sizeof(ECPoint), ctx, pd, cq);
    RDMAConnection conn4 = setup_connection(RECEIVER4_IP, chunk4, (total - 3*quarter) * sizeof(ECPoint), ctx, pd, cq);

    std::cout << "Chunks sent to all receivers. Waiting for results back via RDMA write or TCP...\n";

    // Allocate GPU memory for receiving the results
    ECPoint* recv_buf1 = nullptr;
    ECPoint* recv_buf2 = nullptr;
    ECPoint* recv_buf3 = nullptr;
    ECPoint* recv_buf4 = nullptr;
    
    cudaError_t recv_err1 = cudaMalloc((void**)&recv_buf1, sizeof(ECPoint));
    cudaError_t recv_err2 = cudaMalloc((void**)&recv_buf2, sizeof(ECPoint));
    cudaError_t recv_err3 = cudaMalloc((void**)&recv_buf3, sizeof(ECPoint));
    cudaError_t recv_err4 = cudaMalloc((void**)&recv_buf4, sizeof(ECPoint));
    
    if (recv_err1 != cudaSuccess || recv_err2 != cudaSuccess || 
        recv_err3 != cudaSuccess || recv_err4 != cudaSuccess) {
        std::cerr << "cudaMalloc failed for receive buffers: " 
                  << cudaGetErrorString(recv_err1) << ", "
                  << cudaGetErrorString(recv_err2) << ", "
                  << cudaGetErrorString(recv_err3) << ", "
                  << cudaGetErrorString(recv_err4) << std::endl;
        return -1;
    }
    
    // Register GPU memory with RDMA for all receivers
    ibv_mr *recv_mr1 = ibv_reg_mr(pd, recv_buf1, sizeof(ECPoint), 
                                IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_ZERO_BASED);
    ibv_mr *recv_mr2 = ibv_reg_mr(pd, recv_buf2, sizeof(ECPoint), 
                                IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_ZERO_BASED);
    ibv_mr *recv_mr3 = ibv_reg_mr(pd, recv_buf3, sizeof(ECPoint), 
                                IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_ZERO_BASED);
    ibv_mr *recv_mr4 = ibv_reg_mr(pd, recv_buf4, sizeof(ECPoint), 
                                IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_ZERO_BASED);
    
    if (!recv_mr1 || !recv_mr2 || !recv_mr3 || !recv_mr4) {
        std::cerr << "Failed to register memory regions for receiving results" << std::endl;
        cudaFree(recv_buf1);
        cudaFree(recv_buf2);
        cudaFree(recv_buf3);
        cudaFree(recv_buf4);
        return -1;
    }

    // Send our receive buffer info to each receiver using existing connections
    rdma_info recv_info1 = {
        .addr = (uint64_t)(uintptr_t)recv_buf1,
        .rkey = recv_mr1->rkey
    };
    rdma_info recv_info2 = {
        .addr = (uint64_t)(uintptr_t)recv_buf2,
        .rkey = recv_mr2->rkey
    };
    rdma_info recv_info3 = {
        .addr = (uint64_t)(uintptr_t)recv_buf3,
        .rkey = recv_mr3->rkey
    };
    rdma_info recv_info4 = {
        .addr = (uint64_t)(uintptr_t)recv_buf4,
        .rkey = recv_mr4->rkey
    };

    // Send receive buffer info to receivers using existing connections
    ssize_t send_bytes1 = send(conn1.sockfd, &recv_info1, sizeof(recv_info1), 0);
    ssize_t send_bytes2 = send(conn2.sockfd, &recv_info2, sizeof(recv_info2), 0);
    ssize_t send_bytes3 = send(conn3.sockfd, &recv_info3, sizeof(recv_info3), 0);
    ssize_t send_bytes4 = send(conn4.sockfd, &recv_info4, sizeof(recv_info4), 0);
    
    // Receive QP info from receivers
    rdma_info receiver1_info, receiver2_info, receiver3_info, receiver4_info;
    ssize_t recv_bytes1 = recv(conn1.sockfd, &receiver1_info, sizeof(receiver1_info), 0);
    ssize_t recv_bytes2 = recv(conn2.sockfd, &receiver2_info, sizeof(receiver2_info), 0);
    ssize_t recv_bytes3 = recv(conn3.sockfd, &receiver3_info, sizeof(receiver3_info), 0);
    ssize_t recv_bytes4 = recv(conn4.sockfd, &receiver4_info, sizeof(receiver4_info), 0);

    // Create QPs for receiving results
    ibv_cq *recv_cq1 = ibv_create_cq(ctx, 10, nullptr, nullptr, 0);
    ibv_cq *recv_cq2 = ibv_create_cq(ctx, 10, nullptr, nullptr, 0);
    ibv_cq *recv_cq3 = ibv_create_cq(ctx, 10, nullptr, nullptr, 0);
    ibv_cq *recv_cq4 = ibv_create_cq(ctx, 10, nullptr, nullptr, 0);
    
    if (!recv_cq1 || !recv_cq2 || !recv_cq3 || !recv_cq4) {
        std::cerr << "Failed to create completion queues: " << strerror(errno) << std::endl;
        if (recv_cq1) ibv_destroy_cq(recv_cq1);
        if (recv_cq2) ibv_destroy_cq(recv_cq2);
        if (recv_cq3) ibv_destroy_cq(recv_cq3);
        if (recv_cq4) ibv_destroy_cq(recv_cq4);
        return -1;
    }
    
    ibv_qp_init_attr recv_qp_init_attr = {};
    recv_qp_init_attr.send_cq = recv_cq1;
    recv_qp_init_attr.recv_cq = recv_cq1;
    recv_qp_init_attr.qp_type = IBV_QPT_RC;
    recv_qp_init_attr.cap.max_send_wr = 1;
    recv_qp_init_attr.cap.max_recv_wr = 1;
    recv_qp_init_attr.cap.max_send_sge = 1;
    recv_qp_init_attr.cap.max_recv_sge = 1;
    
    ibv_qp *recv_qp1 = ibv_create_qp(pd, &recv_qp_init_attr);
    recv_qp_init_attr.send_cq = recv_cq2;
    recv_qp_init_attr.recv_cq = recv_cq2;
    ibv_qp *recv_qp2 = ibv_create_qp(pd, &recv_qp_init_attr);
    recv_qp_init_attr.send_cq = recv_cq3;
    recv_qp_init_attr.recv_cq = recv_cq3;
    ibv_qp *recv_qp3 = ibv_create_qp(pd, &recv_qp_init_attr);
    recv_qp_init_attr.send_cq = recv_cq4;
    recv_qp_init_attr.recv_cq = recv_cq4;
    ibv_qp *recv_qp4 = ibv_create_qp(pd, &recv_qp_init_attr);
    
    if (!recv_qp1 || !recv_qp2 || !recv_qp3 || !recv_qp4) {
        std::cerr << "Failed to create QPs: " << strerror(errno) << std::endl;
        if (recv_qp1) ibv_destroy_qp(recv_qp1);
        if (recv_qp2) ibv_destroy_qp(recv_qp2);
        if (recv_qp3) ibv_destroy_qp(recv_qp3);
        if (recv_qp4) ibv_destroy_qp(recv_qp4);
        ibv_destroy_cq(recv_cq1);
        ibv_destroy_cq(recv_cq2);
        ibv_destroy_cq(recv_cq3);
        ibv_destroy_cq(recv_cq4);
        return -1;
    }
    
    // Modify QPs to RTS
    modify_qp_to_rts(recv_qp1, receiver1_info.qp_num, receiver1_info.lid, receiver1_info.gid);
    modify_qp_to_rts(recv_qp2, receiver2_info.qp_num, receiver2_info.lid, receiver2_info.gid);
    modify_qp_to_rts(recv_qp3, receiver3_info.qp_num, receiver3_info.lid, receiver3_info.gid);
    modify_qp_to_rts(recv_qp4, receiver4_info.qp_num, receiver4_info.lid, receiver4_info.gid);

    // Wait for results to be written by receivers
    // In a real implementation, you'd want to poll or use completion events
    sleep(2);  // Simple wait for demonstration

    // Copy results back to host memory for printing
    ECPoint result1, result2, result3, result4;
    cudaError_t copy_err1 = cudaMemcpy(&result1, recv_buf1, sizeof(ECPoint), cudaMemcpyDeviceToHost);
    cudaError_t copy_err2 = cudaMemcpy(&result2, recv_buf2, sizeof(ECPoint), cudaMemcpyDeviceToHost);
    cudaError_t copy_err3 = cudaMemcpy(&result3, recv_buf3, sizeof(ECPoint), cudaMemcpyDeviceToHost);
    cudaError_t copy_err4 = cudaMemcpy(&result4, recv_buf4, sizeof(ECPoint), cudaMemcpyDeviceToHost);
    
    if (copy_err1 != cudaSuccess || copy_err2 != cudaSuccess || 
        copy_err3 != cudaSuccess || copy_err4 != cudaSuccess) {
        std::cerr << "Failed to copy results from device: " 
                  << cudaGetErrorString(copy_err1) << ", "
                  << cudaGetErrorString(copy_err2) << ", "
                  << cudaGetErrorString(copy_err3) << ", "
                  << cudaGetErrorString(copy_err4) << std::endl;
        return -1;
    }

    std::cout << "Result from receiver 1: (" << result1.x << ", " << result1.y << ")\n";
    std::cout << "Result from receiver 2: (" << result2.x << ", " << result2.y << ")\n";
    std::cout << "Result from receiver 3: (" << result3.x << ", " << result3.y << ")\n";
    std::cout << "Result from receiver 4: (" << result4.x << ", " << result4.y << ")\n";

    const int field_p = 11;
    const int curve_a = 3;
    const int windowSize = 3;

    ECPoint sum = {0,0,1};

    sum = ecAddTrueCPU(sum, result1, curve_a, field_p);
    sum = ecAddTrueCPU(sum, result2, curve_a, field_p);
    sum = ecAddTrueCPU(sum, result3, curve_a, field_p);
    sum = ecAddTrueCPU(sum, result4, curve_a, field_p);

    std::cout << "Final point: (" << sum.x << ", " << sum.y << ")\n";

    // Cleanup
    ibv_dereg_mr(recv_mr1);
    ibv_dereg_mr(recv_mr2);
    ibv_dereg_mr(recv_mr3);
    ibv_dereg_mr(recv_mr4);
    ibv_destroy_qp(recv_qp1);
    ibv_destroy_qp(recv_qp2);
    ibv_destroy_qp(recv_qp3);
    ibv_destroy_qp(recv_qp4);
    ibv_destroy_cq(recv_cq1);
    ibv_destroy_cq(recv_cq2);
    ibv_destroy_cq(recv_cq3);
    ibv_destroy_cq(recv_cq4);
    ibv_destroy_qp(conn1.qp);
    ibv_destroy_qp(conn2.qp);
    ibv_destroy_qp(conn3.qp);
    ibv_destroy_qp(conn4.qp);
    ibv_destroy_cq(cq);
    ibv_dealloc_pd(pd);
    ibv_close_device(ctx);
    ibv_free_device_list(dev_list);

    close(conn1.sockfd);
    close(conn2.sockfd);
    close(conn3.sockfd);
    close(conn4.sockfd);

    free(chunk1);
    free(chunk2);
    free(chunk3);
    free(chunk4);
    cudaFree(recv_buf1);
    cudaFree(recv_buf2);
    cudaFree(recv_buf3);
    cudaFree(recv_buf4);

    return 0;
}
