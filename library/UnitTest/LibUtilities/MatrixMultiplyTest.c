#include "LibUtiltiesTest.h"

int main(void){
    int M=3, N=4, K=3;
    double *A = BuildVector(M*N);
    double *B = BuildVector(N*K);
    double *C = BuildVector(M*K);
    double *exctC = BuildVector(M*K);

    double tempA[3][4] = {{3,4,1,2},
                          {2,5,8,3},
                          {7,4,2,4}};

    double tempB[4][3] = {{2,4,-4},
                          {1,0,5},
                          {0,3,3},
                          {-5,2,2}};

    double tempC[3][3] = {{0, 19, 15},
                          {-6, 38, 47},
                          {-2, 42, 6}};
    int i,j,info;

    int sk = 0;
    for(i=0;i<M;i++){
        for(j=0;j<N;j++)
            A[sk++] = tempA[i][j];
    }

    sk = 0;
    for(i=0;i<N;i++){
        for(j=0;j<K;j++)
            B[sk++] = tempB[i][j];
    }

    sk = 0;
    for(i=0;i<M;i++){
        for(j=0;j<K;j++) {
            exctC[sk++] = tempC[i][j];
        }
    }

    dgemm_(M, M, K, N, A, B, C);

    info = CreateVectorTest("Matrix Inverse", C, exctC, M*K);

    DestroyVector(A);
    DestroyVector(B);
    DestroyVector(C);
    DestroyVector(exctC);
    return info;
}
