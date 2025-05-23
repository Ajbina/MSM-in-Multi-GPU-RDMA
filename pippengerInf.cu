#include "pippengerInf.h"
#include <cuda_runtime.h>
#include <cuda.h>
#include <stdio.h>

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

/*
__global__ void pippengerSumming(ECPoint *result, const int blocksNeeded){
    const int field_p = 11;
    const int curve_a = 3;
    ECPoint sum = {0,0,1};
    printf("result: %d\n", blocksNeeded);
    for(int i = 0; i < blocksNeeded; i++){
        sum = ecAddTrue(sum, result[i], curve_a, field_p);
        printf("value adding: %d, %d\n", result[i].x, result[i].y);
    }

    result[0] = sum;
}*/

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

    // if(tid == 0){
    //     printf("TESTING more: %d, %d, %d\n", 1,blockIdx.x, numBucketsTotal);

    // }
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
            //printf("windowSum: x = %d, y = %d\n", windowSum.x, windowSum.y);
        }

        result[blockIdx.x] = windowSum;
        //printf("--------------- block num: %d, num: %d sum: %d\n", blockIdx.x, 0, result[0].x);
        //printf("--------------- block num: %d, num: %d sum: %d\n", blockIdx.x, 1, result[1].x);
    }
    __syncthreads();
}
