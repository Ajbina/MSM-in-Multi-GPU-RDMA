#include <iostream>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <chrono>
#include <random>
#include <string>
#include <sstream>
#include <fstream>
  
#include "pippengerInf.h"  

const int field_p = 11;
const int curve_a = 3;
const int windowSize = 3;

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

ECPoint MSMBadApproach(std::vector<ECPoint> pointsIn, std::vector<int> scalarsIn){
    ECPoint finalPoint = {-1,-1,1};
    ECPoint sum = {-1,-1,1};
    int pointAmount = pointsIn.size();
    for(int i = 0; i < pointAmount; i++){
        int digit = scalarsIn[i];
        sum.x = pointsIn[i].x;
        sum.y = pointsIn[i].y;
        sum.infinityFlag = pointsIn[i].infinityFlag;

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
            sum = ecDoubleTrueCPU(sum, curve_a, field_p);
        }
        for(int i = 0; i < remaining; i++){
            sum = ecAddTrueCPU(sum, singleVal, curve_a, field_p);
        }

        //printf("point: (%d,%d)\n",sum.x,sum.y);
        finalPoint = ecAddTrueCPU(finalPoint, sum, curve_a, field_p);
    }
    return finalPoint;
}


int main() {

    //std::vector<int> h_scalars = {1, 12341234, 3, 4, 5, 6, 232134, 7, 7, 1001, 1234};

    /*std::vector<ECPoint> h_points = {
        {3, 4, 0},
        {3, 7, 0},
        {1, 10, 0},
        {5, 9, 0},
        {5, 2, 0},
        {1, 1, 0},
        {4, 8, 0},
        {4, 3, 0},
        {9, 6, 0},
        {9, 5, 0},
        {0, 0, 0}
    };*/

    /*
    std::vector<ECPoint> h_points = {
        {0, 0, 0}, {1, 2, 0}, {1, 9, 0}, {2, 5, 0}, {2, 6, 0}, {3, 5, 0}, {3, 6, 0}, {6, 5, 0}, {6, 6, 0}, {7, 1, 0},
        {7, 10, 0}, {0, 0, 0}, {1, 2, 0}, {1, 9, 0}, {2, 5, 0}, {2, 6, 0}, {3, 5, 0}, {3, 6, 0}, {6, 5, 0}, {6, 6, 0},
        {7, 1, 0}, {7, 10, 0}, {0, 0, 0}, {1, 2, 0}, {1, 9, 0}, {2, 5, 0}, {2, 6, 0}, {3, 5, 0}, {3, 6, 0}, {6, 5, 0},
        {6, 6, 0}, {7, 1, 0}, {7, 10, 0}, {0, 0, 0}, {1, 2, 0}, {1, 9, 0}, {2, 5, 0}, {2, 6, 0}, {3, 5, 0}, {3, 6, 0},
        {6, 5, 0}, {6, 6, 0}, {7, 1, 0}, {7, 10, 0}, {0, 0, 0}, {1, 2, 0}, {1, 9, 0}, {2, 5, 0}, {2, 6, 0}, {3, 5, 0},
        {3, 6, 0}, {6, 5, 0}, {6, 6, 0}, {7, 1, 0}, {7, 10, 0}, {0, 0, 0}, {1, 2, 0}, {1, 9, 0}, {2, 5, 0}, {2, 6, 0},
        {3, 5, 0}, {3, 6, 0}, {6, 5, 0}, {6, 6, 0}, {7, 1, 0}, {7, 10, 0}, {0, 0, 0}, {1, 2, 0}, {1, 9, 0}, {2, 5, 0},
        {2, 6, 0}, {3, 5, 0}, {3, 6, 0}, {6, 5, 0}, {6, 6, 0}, {7, 1, 0}, {7, 10, 0}, {0, 0, 0}, {1, 2, 0}, {1, 9, 0},
        {2, 5, 0}, {2, 6, 0}, {3, 5, 0}, {3, 6, 0}, {6, 5, 0}, {6, 6, 0}, {7, 1, 0}, {7, 10, 0}, {0, 0, 0}, {1, 2, 0},
        {1, 9, 0}, {2, 5, 0}, {2, 6, 0}, {3, 5, 0}, {3, 6, 0}, {6, 5, 0}, {6, 6, 0}, {7, 1, 0}, {7, 10, 0}, {0, 0, 0},
        {1, 2, 0}, {1, 9, 0}, {2, 5, 0}, {2, 6, 0}, {3, 5, 0}, {3, 6, 0}, {6, 5, 0}, {6, 6, 0}, {7, 1, 0}, {7, 10, 0},
        {0, 0, 0}, {1, 2, 0}, {1, 9, 0}, {2, 5, 0}, {2, 6, 0}, {3, 5, 0}, {3, 6, 0}, {6, 5, 0}, {6, 6, 0}, {7, 1, 0},
        {7, 10, 0}, {0, 0, 0}, {1, 2, 0}, {1, 9, 0}, {2, 5, 0}, {2, 6, 0}, {3, 5, 0}, {3, 6, 0}, {6, 5, 0}, {6, 6, 0},
        {7, 1, 0}, {7, 10, 0}, {0, 0, 0}, {1, 2, 0}, {1, 9, 0}, {2, 5, 0}, {2, 6, 0}, {3, 5, 0}, {3, 6, 0}, {6, 5, 0},
        {6, 6, 0}, {7, 1, 0}, {7, 10, 0}, {0, 0, 0}, {1, 2, 0}, {1, 9, 0}, {2, 5, 0}, {2, 6, 0}, {3, 5, 0}, {3, 6, 0},
        {6, 5, 0}, {6, 6, 0}, {7, 1, 0}, {7, 10, 0}, {0, 0, 0}, {1, 2, 0}, {1, 9, 0}, {2, 5, 0}, {2, 6, 0}, {3, 5, 0},
        {3, 6, 0}, {6, 5, 0}, {6, 6, 0}, {7, 1, 0}, {7, 10, 0}, {0, 0, 0}, {1, 2, 0}, {1, 9, 0}, {2, 5, 0}, {2, 6, 0},
        {3, 5, 0}, {3, 6, 0}, {6, 5, 0}, {6, 6, 0}, {7, 1, 0}, {7, 10, 0}, {0, 0, 0}, {1, 2, 0}, {1, 9, 0}, {2, 5, 0},
        {2, 6, 0}, {3, 5, 0}, {3, 6, 0}, {6, 5, 0}, {6, 6, 0}, {7, 1, 0}, {7, 10, 0}, {0, 0, 0}, {1, 2, 0}, {1, 9, 0},
        {2, 5, 0}, {2, 6, 0}, {3, 5, 0}, {3, 6, 0}, {6, 5, 0}, {6, 6, 0}, {7, 1, 0}, {7, 10, 0}, {0, 0, 0}, {1, 2, 0},
        {1, 9, 0}, {2, 5, 0}, {2, 6, 0}, {3, 5, 0}, {3, 6, 0}, {6, 5, 0}, {6, 6, 0}, {7, 1, 0}, {7, 10, 0}, {0, 0, 0},
        {1, 2, 0}, {1, 9, 0}, {2, 5, 0}, {2, 6, 0}, {3, 5, 0}, {3, 6, 0}, {6, 5, 0}, {6, 6, 0}, {7, 1, 0}, {7, 10, 0}  
    };*/
        // 0,0
        // 1,1
        // 10,1
        // 4,3
        // 7,3
        // 8,4
        // 3,4
        // 2,5
        // 9,5
        // 6,9
        // 5,9

    //h_scalars = {1, 2, 3, 4, 5, 6, 7, 7, 7, 11};//, 11};

    std::string filename = "data/points.csv";
    std::ifstream file(filename);

    std::vector<ECPoint> h_points;

    std::string line;
    std::getline(file, line);
    int numPointsCSV = std::stoi(line);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<int> coords;

        while (std::getline(ss, token, ',')) {
            coords.push_back(std::stoi(token));
        }
        ECPoint point;
        point.x = coords[0];
        point.y = coords[1];
        point.infinityFlag = coords[2];

        h_points.emplace_back(point);
    }

    file.close();

    const int numScalars = numPointsCSV;
    const int lowerBound = 1;
    const int upperBound = 1000000;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(lowerBound, upperBound);

    std::vector<int> h_scalars(numScalars);
    for (int& scalar : h_scalars) {
        scalar = dis(gen);
    }
    /*
    std::vector<ECPoint> h_points = {
        {0, 0, 0},
        {1, 2, 0},
        {1, 9, 0},
        {2, 5, 0},
        {2, 6, 0},
        {3, 5, 0},
        {3, 6, 0},
        {6, 5, 0},
        {6, 6, 0},
        {7, 1, 0},
        {7, 10, 0}
    };*/

    if (h_scalars.size() != h_points.size()) {
        printf("Scalars Total: %d, Points total: %d\n", h_scalars.size(), h_points.size());
        std::cerr << "Error: The number of scalars and EC points aint it" << std::endl;
        return -1;
    }

    const int numElements = h_points.size();
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
    ECPoint *d_points = nullptr;     
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
    //    blocksNeeded = 1;
    //}
    
    cudaMalloc(&d_scalars, numElements * sizeof(int));
    cudaMalloc(&d_points, numElements * sizeof(ECPoint));
    cudaMalloc(&d_result, blocksNeeded * sizeof(ECPoint));
    //cudaMalloc(&finalResult, 1 * sizeof(ECPoint));
    cudaMalloc(&allLocalBuckets, numBucketsTotal * numElements * sizeof(ECPoint));//maybe numBucketsTotal * sizeof(ECPoint) * numElements 

    cudaMemcpy(d_scalars, h_scalars.data(), numElements * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_points, h_points.data(), numElements * sizeof(ECPoint), cudaMemcpyHostToDevice);

    //size_t sharedMemSize = numBucketsTotal * sizeof(ECPoint) * numElements + numBucketsTotal * sizeof(ECPoint) + 5000+5000; // windowsize = 3 - (1<<3) = 7
    //size_t sharedMemSize = numBucketsTotal * sizeof(ECPoint);
    size_t sharedMemSize = 80;
    
    /*int threadsNeeded = 0;
    if(numBucketsTotal > numElements){
        threadsNeeded = numBucketsTotal;
    } else {
        threadsNeeded = numElements;
    }*/
    std::cout << "start" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    pippengerKernel<<<blocksNeeded, 1024, sharedMemSize>>>(d_scalars, numElements, d_points, d_result, numBucketsTotal, allLocalBuckets);
    
    cudaDeviceSynchronize();

    //pippengerSumming(const int numElements, ECPoint *result, const int numBucketsTotal, ECPoint *allLocalBuckets)

    //pippengerSumming<<<1, 1024>>>(numElements, d_result, numBucketsTotal, allLocalBuckets);

    //cudaDeviceSynchronize();



    std::vector<ECPoint> h_result(blocksNeeded);
    cudaMemcpy(h_result.data(), d_result, sizeof(ECPoint), cudaMemcpyDeviceToHost);
    
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Final result: (" << h_result[0].x << ", " << h_result[0].y << ")" << std::endl;

    cudaFree(d_scalars);
    cudaFree(d_points);
    cudaFree(d_result);
    cudaFree(allLocalBuckets);

    //if (pointInVector(h_points, h_result[0].x, h_result[0].y)) {
    //    std::cout << "Wow u did it" << std::endl;
    //} else {
    //    std::cout << "you have failed" << std::endl;
    //}
    //std::chrono::duration<double> elapsed_seconds = end - start;
    //std::cout << "Time GPU: " << elapsed_seconds.count() << " seconds" << std::endl;
    //slow msm approach
    
    auto start2 = std::chrono::high_resolution_clock::now();
    ECPoint finalPoint = MSMBadApproach(h_points, h_scalars);
    auto end2 = std::chrono::high_resolution_clock::now();
    
    std::cout << "Final result: (" << finalPoint.x << ", " << finalPoint.y << ")" << std::endl;


    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Time GPU: " << elapsed_seconds.count() << " seconds" << std::endl;
    std::chrono::duration<double> elapsed_seconds2 = end2 - start2;
    std::cout << "Time CPU: " << elapsed_seconds2.count() << " seconds" << std::endl;
    
    return 0;
}