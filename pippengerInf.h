#ifndef PIPPENGERINF_H
#define PIPPENGERINF_H

#ifdef __cplusplus
extern "C"{
    #endif

    //struct for an elliptic curve point
    struct ECPoint{
        int x;
        int y;
        int infinityFlag;
        //ECPoint(int x_, int y_, int infinityFlag_) : x(x_), y(y_), infinityFlag(infinityFlag_) {}
    };

    //CUDA Kernel declaration for puppenger bucket
    //takes as input array of scalars/ec points and outputs results
    __global__ void pippengerKernel(const int *scalars, const int numElements, const ECPoint *points, ECPoint *result, int numBucketsTotal, ECPoint *allLocalBuckets);
    __global__ void pippengerSumming(const int numElements, ECPoint *result, const int numBucketsTotal, ECPoint *allLocalBuckets);

    //Device function prototypes for elliptic curve arithmetic
    __device__ ECPoint ecAddTrue(const ECPoint &a, const ECPoint &b, int a_coeff, int p);
    __device__ ECPoint ecDoubleTrue(const ECPoint &p, int a_coeff, int field_p);    

    #ifdef __cplusplus
}
#endif

#endif