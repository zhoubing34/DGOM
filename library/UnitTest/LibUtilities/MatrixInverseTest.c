#include "LibUtiltiesTest.h"

/**
 * @brief
 * Test of Matrix inverse function
 */
int main(void){
    int i,j,N=3,sk=0, info;
    double Ct[3][3] = {{0.25, 0.5, 0.5},
                      {0, 1.0, 0.25},
                      {0.25, 0.5, 0.25}};
    double invCt[3][3] = {{-2, -2, 6},
                         {-1, 1, 1},
                         {4, 0, -4}};
    double **C    = Matrix_create(N, N);
    double **invC = Matrix_create(N, N);
    double *temp  = Vector_create(N * N);



    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            C[i][j] = Ct[i][j];
            invC[i][j] = invCt[i][j];
            temp[sk++] = C[i][j];
        }
    }

    invM(temp, N);

    sk = 0;
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
             C[i][j] = temp[sk++];
        }
    }

    info = CreateMatrixTest("Matrix Inverse", C, invC, N, N);

    Matrix_free(C);
    Matrix_free(invC);
    Vector_free(temp);

    return info;
}
