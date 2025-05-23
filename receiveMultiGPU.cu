#include <iostream>
#include <chrono>
#include <random>
#include <cuda_runtime.h>
#include <iomanip>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <unistd.h>
#include <arpa/inet.h>
#include <infiniband/verbs.h>
#include <string>
#include <sstream>
#include <fstream>

//ur going to need to do the data stuff in the send to make this work and adjust the receive amount

#define TCP_PORT 18515
#define SENDER_IP "128.105.144.57" // Dummy IP for now, will be replaced with actual IP

struct rdma_info {
 uint64_t addr;
 uint32_t rkey;
 uint32_t qp_num;
 uint16_t lid;
 uint8_t gid[16];
};

struct ECPoint{
 int x;
 int y;
 int infinityFlag;
};

//const int field_p = 11;
//const int curve_a = 3;
const int windowSize = 3;

__device__ int mod(int a, int m) {
 int r = a % m;
 return (r < 0) ? r + m : r;
} 

__device__ int extractWindow(int scalar, int windowIndex, int windowSize){
 int mask = (1 << windowSize) - 1; //if windowSize == 3, mask = 0x111 = 7
 return (scalar >> (windowIndex * windowSize)) & mask;
}

//device function for inverse
__device__ int modInverse(int a, int p) {
 a = mod(a, p);
 for (int x = 1; x < p; x++) {
 if (mod(a * x, p) == 1) {
 return x;
 }
 }
 return 0;
}

//elliptic curve point doubling
__device__ ECPoint ecDoubleTrue(const ECPoint &P, int a, int p){
 ECPoint R;
 //Special case: if P.y == 0, then doubling results in point at infinity
 if (P.y == 0){
 return {-1, -1, 1};
 }
 int infP = *((volatile int*)&(P.infinityFlag));
 if (infP == 1){
 return {-1, -1, 1};
 }

 int numerator = mod(3 * P.x * P.x + a, p);
 int denominator = mod(2 * P.y, p);
 int invDenom = modInverse(denominator, p);
 int lambda = mod(numerator * invDenom, p);

 R.x = mod(lambda * lambda - 2 * P.x, p);
 R.y = mod(lambda * (P.x - R.x) - P.y, p);
 R.infinityFlag = 0;
 return R;
}

//elliptic curve point addition
__device__ ECPoint ecAddTrue(const ECPoint &P, const ECPoint &Q, int a, int p) {
 int infP = *((volatile int*)&(P.infinityFlag));
 int infQ = *((volatile int*)&(Q.infinityFlag));
 ECPoint R;
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
 R = ecDoubleTrue(P, a, p);
 return R;
 }

 int numerator = mod(Q.y - P.y, p);
 int denominator = mod(Q.x - P.x, p);
 int invDenom = modInverse(denominator, p);
 int lambda = mod(numerator * invDenom, p);

 R.x = mod(lambda * lambda - P.x - Q.x, p);
 R.y = mod(lambda * (P.x - R.x) - P.y, p);
 R.infinityFlag = 0;
 return R;
}

__global__ void pippengerKernel(const int *scalars, const int numElements, const ECPoint *points, ECPoint *result, const int numBucketsTotal, ECPoint *allLocalBuckets){
 const int field_p = 11;
 const int curve_a = 3;
 const int windowSize = 3;
 const int bucketsPerWindow = 1 << windowSize;
 const int numWindows = numBucketsTotal / bucketsPerWindow;
 
 int tid = threadIdx.x;

 int idx = blockIdx.x * blockDim.x + threadIdx.x;
 
 if (tid == 0) {
 result[blockIdx.x].x = 0;
 result[blockIdx.x].y = 0;
 result[blockIdx.x].infinityFlag = 0;
 }
 __syncthreads();
 
 ECPoint localBuckets[80];//80 covers upt to 2^30 scalars with window size 3 cuz there r 10 windows of buckets of 8
 for (int b = 0; b < numBucketsTotal; b++){
 localBuckets[b].x = 0;
 localBuckets[b].y = 0;
 localBuckets[b].infinityFlag = 1;
 }
 int mask = (1 << windowSize) - 1;
 
 if (idx < numElements){
 int scalar = int(scalars[idx]);

 for (int w = 0; w < numWindows; w++){
 int shiftAmount = w * windowSize;
 int digit = (scalar >> shiftAmount) & mask;
 if (digit != 0) {
 int bucketIdx = digit;
 localBuckets[bucketIdx + bucketsPerWindow * w] = ecAddTrue(localBuckets[bucketIdx + bucketsPerWindow * w], points[idx], curve_a, field_p);
 }
 }
 }

 if(idx < numElements){
 for (int b = 1; b < numBucketsTotal; b++){
 allLocalBuckets[idx * numBucketsTotal + b] = localBuckets[b];
 }
 }
 __shared__ ECPoint buckets[80];
 
 if(tid == 0){
 for(int i = 0; i < 80; i++){
 buckets[i].x = 0;
 buckets[i].y = 0;
 buckets[i].infinityFlag = 1;
 }
 }
 __syncthreads();

 if(tid < numBucketsTotal){
 int b = tid;
 int digit = b % bucketsPerWindow;
 ECPoint sum = {0,0,1};
 if(digit == 0){
 buckets[b] = sum;
 } else {
 for (int t = 0; t < numElements; t++){
 sum = ecAddTrue(sum, allLocalBuckets[t * numBucketsTotal + tid], curve_a, field_p);
 }
 int inf = *((volatile int*)&(sum.infinityFlag));
 if(inf == 0 && digit != 0){
 int n = digit;
 int counter = 0;
 int current = 2;

 while (current <= n) {
 counter++;
 current *= 2;
 }

 current /= 2;
 int remaining = n - current;

 ECPoint singleVal = {sum.x,sum.y};
 for(int i = 0; i < counter; i++){
 sum = ecDoubleTrue(sum, curve_a, field_p);
 }
 for(int i = 0; i < remaining; i++){
 sum = ecAddTrue(sum, singleVal, curve_a, field_p);
 }
 }
 buckets[b] = sum;
 }
 }
 __syncthreads();

 if(tid < numBucketsTotal){
 ECPoint windowSum = buckets[tid];

 int windowScaleMultiplier = tid / bucketsPerWindow;
 for (int i = 0; i < windowScaleMultiplier * windowSize; i++){
 windowSum = ecDoubleTrue(windowSum, curve_a, field_p);
 }

 buckets[tid] = windowSum;
 }

 if (tid == 0) {

 ECPoint windowSum = {0,0,1};

 for(int i = 0; i < numBucketsTotal; i++){
 windowSum = ecAddTrue(buckets[i], windowSum, curve_a, field_p);
 }

 result[blockIdx.x] = windowSum;
 }
 __syncthreads();
}

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

int main() {
 std::cout << "Receiver: RDMA setup starting...\n";

 ibv_device **dev_list = ibv_get_device_list(nullptr);
 ibv_context *ctx = ibv_open_device(dev_list[0]);
 ibv_pd *pd = ibv_alloc_pd(ctx);

 //int *recv_buf = static_cast<int*>(malloc(2 * sizeof(int)));
 //memset(recv_buf, 0, 2 * sizeof(int));
 ECPoint* recv_buf = nullptr;
 //cudaMalloc((void**)&recv_buf, 2 * sizeof(int));




// CHANGE HERE
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 int numberAmountRec = 8388608;
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




 size_t bytesRec = numberAmountRec * sizeof(ECPoint);
 cudaError_t err = cudaMalloc((void **)&recv_buf, bytesRec);
 if (err != cudaSuccess) {
 fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(err));
 exit(1);
 }

 ibv_mr *mr = ibv_reg_mr(pd, recv_buf, bytesRec,
 IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_ZERO_BASED);

 ibv_cq *cq = ibv_create_cq(ctx, 10, nullptr, nullptr, 0);

 ibv_qp_init_attr qp_init_attr = {};
 qp_init_attr.send_cq = cq;
 qp_init_attr.recv_cq = cq;
 qp_init_attr.qp_type = IBV_QPT_RC;
 qp_init_attr.cap.max_send_wr = 1;
 qp_init_attr.cap.max_recv_wr = 1;
 qp_init_attr.cap.max_send_sge = 1;
 qp_init_attr.cap.max_recv_sge = 1;

 ibv_qp *qp = ibv_create_qp(pd, &qp_init_attr);

 // Get local LID and GID
 ibv_port_attr port_attr;
 ibv_query_port(ctx, 1, &port_attr);
 uint16_t local_lid = port_attr.lid;
 
 union ibv_gid gid;
 ibv_query_gid(ctx, 1, 0, &gid);

 std::cout << "Receiver: RDMA setup complete. Waiting for metadata...\n";

 // TCP socket to receive metadata
 int sockfd = socket(AF_INET, SOCK_STREAM, 0);
 sockaddr_in servaddr{};
 servaddr.sin_family = AF_INET;
 servaddr.sin_port = htons(TCP_PORT);
 servaddr.sin_addr.s_addr = INADDR_ANY;

 bind(sockfd, (sockaddr*)&servaddr, sizeof(servaddr));
 listen(sockfd, 1);
 int connfd = accept(sockfd, nullptr, nullptr);

 // First receive sender's QP info
 rdma_info remote_info;
 recv(connfd, &remote_info, sizeof(remote_info), 0);
 
 // Then send our QP info
 rdma_info local_info = {
 .addr = reinterpret_cast<uint64_t>(recv_buf),
 .rkey = mr->rkey,
 .qp_num = qp->qp_num,
 .lid = local_lid
 };
 memcpy(local_info.gid, gid.raw, 16);
 
 send(connfd, &local_info, sizeof(local_info), 0);
// close(connfd);
// close(sockfd);

 // Modify QP to RTS
 modify_qp_to_rts(qp, remote_info.qp_num, remote_info.lid, remote_info.gid);

 std::cout << "Receiver: Metadata exchanged. Remote addr = " << std::hex << remote_info.addr 
 << ", rkey = " << remote_info.rkey << "\n";

 // RDMA read
 ibv_sge sge{};
 std::cout << "testing loc: 1" << std::endl;
 sge.addr = reinterpret_cast<uint64_t>(recv_buf);
 sge.length = bytesRec;
 sge.lkey = mr->lkey;
 std::cout << "testing loc: 2" << std::endl;
 ibv_send_wr wr{}, *bad_wr = nullptr;
 wr.wr_id = 1;
 wr.sg_list = &sge;
 wr.num_sge = 1;
 wr.opcode = IBV_WR_RDMA_READ;
 wr.send_flags = IBV_SEND_SIGNALED;
 wr.wr.rdma.remote_addr = remote_info.addr;
 wr.wr.rdma.rkey = remote_info.rkey;

 //cudaDeviceSynchronize();
 ibv_post_send(qp, &wr, &bad_wr);
 //cudaDeviceSynchronize();
 // Wait for completion
 std::cout << "testing loc: 3" << std::endl;
 ibv_wc wc{};
 while (ibv_poll_cq(cq, 1, &wc) == 0) {}
 
 if (wc.status != IBV_WC_SUCCESS) {
 std::cerr << "RDMA read failed with status: " << ibv_wc_status_str(wc.status) << "\n";
 return 1;
 }


 //std::vector<int> h_scalars = {7, 4, 5, 5, 5, 4, 2, 5, 5, 1};
 

 const int numScalars = numberAmountRec;
 const int lowerBound = 1;
 const int upperBound = 1048576;

 std::random_device rd;
 std::mt19937 gen(rd());
 std::uniform_int_distribution<> dis(lowerBound, upperBound);

 std::vector<int> h_scalars(numScalars);
 for (int& scalar : h_scalars) {
 scalar = dis(gen);
 }



 auto start = std::chrono::high_resolution_clock::now();
 const int numElements = h_scalars.size();
 int max = 0;
 for(int scalarVal : h_scalars){
 if(scalarVal > max){
 max = scalarVal;
 }
 }
 int numZeros = 0;
 for(int i = 31; i >= 0; i--){
 int val = ((max >> i) & 0x0001);
 if(val == 0){
 numZeros++;
 } else {
 break;
 }
 }
 int bitsUsed = 32 - numZeros;
 int numWindows = 0;
 while(numWindows * windowSize < bitsUsed){
 numWindows += 1;
 }

 const int numBucketsTotal = numWindows * (1 << windowSize);

 int *d_scalars = nullptr; 
 //ECPoint *d_points = nullptr; 
 ECPoint *d_result = nullptr; 
 ECPoint *allLocalBuckets = nullptr;
 //ECPoint *finalResult = nullptr;
 int blocksNeeded = 0;
 if(numElements % 1024 != 0){
 blocksNeeded = numElements/1024+1;
 } else {
 blocksNeeded = numElements/1024;
 }
 //if(blocksNeeded == 0){
 // blocksNeeded = 1;
 //}
 
 cudaMalloc(&d_scalars, numElements * sizeof(int));
 //cudaMalloc(&d_points, numElements * sizeof(ECPoint));
 cudaMalloc(&d_result, blocksNeeded * sizeof(ECPoint));
 //cudaMalloc(&finalResult, 1 * sizeof(ECPoint));
 cudaMalloc(&allLocalBuckets, numBucketsTotal * numElements * sizeof(ECPoint));//maybe numBucketsTotal * sizeof(ECPoint) * numElements 

 cudaMemcpy(d_scalars, h_scalars.data(), numElements * sizeof(int), cudaMemcpyHostToDevice);
 //cudaMemcpy(d_points, h_points.data(), numElements * sizeof(ECPoint), cudaMemcpyHostToDevice);

 //size_t sharedMemSize = numBucketsTotal * sizeof(ECPoint) * numElements + numBucketsTotal * sizeof(ECPoint) + 5000+5000; // windowsize = 3 - (1<<3) = 7
 //size_t sharedMemSize = numBucketsTotal * sizeof(ECPoint);
 size_t sharedMemSize = 80;
 
 /*int threadsNeeded = 0;
 if(numBucketsTotal > numElements){
 threadsNeeded = numBucketsTotal;
 } else {
 threadsNeeded = numElements;
 }*/
 //std::cout << "start" << std::endl;
 //auto start = std::chrono::high_resolution_clock::now();
 pippengerKernel<<<blocksNeeded, 1024, sharedMemSize>>>(d_scalars, numElements, recv_buf, d_result, numBucketsTotal, allLocalBuckets);
 
 cudaDeviceSynchronize();
 auto end = std::chrono::high_resolution_clock::now();
 std::vector<ECPoint> h_result(1);
 cudaMemcpy(h_result.data(), d_result, sizeof(ECPoint), cudaMemcpyDeviceToHost);

 std::cout << "Final result: (" << h_result[0].x << ", " << h_result[0].y << ")" << std::endl;

 std::chrono::duration<double> elapsed_seconds = end - start;
 std::cout << "Time GPU: " << elapsed_seconds.count() << " seconds" << std::endl;
 
 // ===== SENDER PART - SENDING RESULT BACK =====
 std::cout << "\nReceiver: Setting up to send result back...\n";
 
 // Receive sender's receive buffer info using existing connection
 rdma_info sender_recv_info;
 ssize_t recv_bytes = recv(connfd, &sender_recv_info, sizeof(sender_recv_info), 0);
 if (recv_bytes != sizeof(sender_recv_info)) {
 std::cerr << "Failed to receive sender's receive buffer info" << std::endl;
 return -1;
 }
 
 // Allocate GPU memory for the result
 ECPoint* send_buf = nullptr;
 cudaError_t send_err = cudaMalloc((void**)&send_buf, sizeof(ECPoint));
 if (send_err != cudaSuccess) {
 std::cerr << "cudaMalloc failed for send buffer: " << cudaGetErrorString(send_err) << std::endl;
 return -1;
 }
 
 // Copy the result to GPU memory
 cudaError_t copy_err = cudaMemcpy(send_buf, &h_result[0], sizeof(ECPoint), cudaMemcpyHostToDevice);
 if (copy_err != cudaSuccess) {
 std::cerr << "cudaMemcpy failed: " << cudaGetErrorString(copy_err) << std::endl;
 cudaFree(send_buf);
 return -1;
 }
 
 // Register GPU memory with RDMA
 ibv_mr *send_mr = ibv_reg_mr(pd, send_buf, sizeof(ECPoint), 
 IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_ZERO_BASED);
 
 if (!send_mr) {
 std::cerr << "Failed to register memory region for sending result: " << strerror(errno) << std::endl;
 cudaFree(send_buf);
 return -1;
 }
 
 // Create a new QP for sending
 ibv_cq *send_cq = ibv_create_cq(ctx, 10, nullptr, nullptr, 0);
 if (!send_cq) {
 std::cerr << "Failed to create completion queue: " << strerror(errno) << std::endl;
 ibv_dereg_mr(send_mr);
 cudaFree(send_buf);
 return -1;
 }
 
 ibv_qp_init_attr send_qp_init_attr = {};
 send_qp_init_attr.send_cq = send_cq;
 send_qp_init_attr.recv_cq = send_cq;
 send_qp_init_attr.qp_type = IBV_QPT_RC;
 send_qp_init_attr.cap.max_send_wr = 1;
 send_qp_init_attr.cap.max_recv_wr = 1;
 send_qp_init_attr.cap.max_send_sge = 1;
 send_qp_init_attr.cap.max_recv_sge = 1;
 
 ibv_qp *send_qp = ibv_create_qp(pd, &send_qp_init_attr);
 if (!send_qp) {
 std::cerr << "Failed to create QP: " << strerror(errno) << std::endl;
 ibv_destroy_cq(send_cq);
 ibv_dereg_mr(send_mr);
 cudaFree(send_buf);
 return -1;
 }
 
 // Send our QP info using existing connection
 rdma_info send_local_info = {
 .addr = (uint64_t)(uintptr_t)send_buf,
 .rkey = send_mr->rkey,
 .qp_num = send_qp->qp_num,
 .lid = local_lid
 };
 memcpy(send_local_info.gid, gid.raw, 16);
 
 ssize_t send_bytes = send(connfd, &send_local_info, sizeof(send_local_info), 0);
 if (send_bytes != sizeof(send_local_info)) {
 std::cerr << "Failed to send QP info: " << strerror(errno) << std::endl;
 ibv_destroy_qp(send_qp);
 ibv_destroy_cq(send_cq);
 ibv_dereg_mr(send_mr);
 cudaFree(send_buf);
 return -1;
 }
 
 // Modify QP to RTS
 modify_qp_to_rts(send_qp, sender_recv_info.qp_num, sender_recv_info.lid, sender_recv_info.gid);
  
 std::cout << "Receiver: Metadata exchanged. Sending result back...\n";
 
 // RDMA write to send the result
 ibv_sge send_sge{};
 send_sge.addr = (uint64_t)(uintptr_t)send_buf;
 send_sge.length = sizeof(ECPoint);
 send_sge.lkey = send_mr->lkey;
 
 ibv_send_wr send_wr{}, *send_bad_wr = nullptr;
 send_wr.wr_id = 2;
 send_wr.sg_list = &send_sge;
 send_wr.num_sge = 1;
 send_wr.opcode = IBV_WR_RDMA_WRITE;
 send_wr.send_flags = IBV_SEND_SIGNALED;
 send_wr.wr.rdma.remote_addr = sender_recv_info.addr;
 send_wr.wr.rdma.rkey = sender_recv_info.rkey;
 
 if (ibv_post_send(send_qp, &send_wr, &send_bad_wr) != 0) {
 std::cerr << "Failed to post send: " << strerror(errno) << std::endl;
 ibv_destroy_qp(send_qp);
 ibv_destroy_cq(send_cq);
 ibv_dereg_mr(send_mr);
 cudaFree(send_buf);
 return -1;
 }
 
 // Wait for completion
 ibv_wc send_wc{};
 int poll_result;
 do {
 poll_result = ibv_poll_cq(send_cq, 1, &send_wc);
 } while (poll_result == 0);
 
 if (poll_result < 0) {
 std::cerr << "Failed to poll CQ: " << strerror(errno) << std::endl;
 ibv_destroy_qp(send_qp);
 ibv_destroy_cq(send_cq);
 ibv_dereg_mr(send_mr);
 cudaFree(send_buf);
 return -1;
 }
 
 if (send_wc.status != IBV_WC_SUCCESS) {
 std::cerr << "RDMA write failed with status: " << ibv_wc_status_str(send_wc.status) << "\n";
 ibv_destroy_qp(send_qp);
 ibv_destroy_cq(send_cq);
 ibv_dereg_mr(send_mr);
 cudaFree(send_buf);
 return -1;
 }
 
 std::cout << "Receiver: Result sent successfully!\n";
 
 // Cleanup
 ibv_dereg_mr(send_mr);
 ibv_destroy_qp(send_qp);
 ibv_destroy_cq(send_cq);
 cudaFree(send_buf);
 
 // Cleanup original receiver resources
 ibv_dereg_mr(mr);
 ibv_destroy_qp(qp);
 ibv_destroy_cq(cq);
 ibv_dealloc_pd(pd);
 ibv_close_device(ctx);
 ibv_free_device_list(dev_list);
 
 close(connfd);
 close(sockfd);
 
 cudaFree(recv_buf);
 cudaFree(d_scalars);
 cudaFree(d_result);
 cudaFree(allLocalBuckets);
 // close(connfd);
// close(sockfd);
 return 0;
} 
