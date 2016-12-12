//
// Created by li12242 on 12/10/16.
//

#ifndef DGOM_MATRIXMULTIPLY_H
#define DGOM_MATRIXMULTIPLY_H

#define M 3
#define K 4
#define N 3

double A[M*K] = {3,4,1,2,
                 2,5,8,3,
                 7,4,2,4};

double B[K*N] = {2,4,-4,
                 1,0,5,
                 0,3,3,
                 -5,2,2};

double C[M*N] = {0, 19, 15,
                 -6, 38, 47,
                 -2, 42, 6};

#endif //DGOM_MATRIXMULTIPLY_H
