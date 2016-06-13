//
// Created by li12242 on 16/6/12.
//

#include <stdio.h>
#include "LibUtilities/LibUtilities.h"
#include "StdRegions/StdRegions.h"

#define Deg 3

int main(int argc, char **argv) {

    StdRegions2d *tri;
    tri = GenStdTriEle(Deg);

    PrintMatrix("Mass matrix", tri->M, tri->Np, tri->Np);
    FreeStdRegions2d(tri);

    return 0;
}
void MatrixMuti(void){
    int i,j,sk=0;
    const int M=Deg;
    double *A, *B, *C;
    A = (doublereal*) malloc(sizeof(double)*Deg*Deg);
    B = (doublereal*) malloc(sizeof(double)*Deg*Deg);
    C = (doublereal*) malloc(sizeof(double)*Deg*Deg);

    sk = 0;
    for(i=0;i<Deg;i++){
        for(j=0;j<Deg;j++){
            A[sk] = i*i + j;
            B[sk] = j*i;
            C[sk] = 0.0;
            sk++;
        }
    }

    sk=0;
    printf("\nA=\n");
    for(i=0;i<Deg;i++){
        for(j=0;j<Deg;j++){
            printf("%lg,", A[sk++]);
        }
        printf("\n");
    }

    sk=0;
    printf("\nB=\n");
    for(i=0;i<Deg;i++){
        for(j=0;j<Deg;j++){
            printf("%lg,", B[sk++]);
        }
        printf("\n");
    }


    dgemm_(Deg, Deg, Deg, Deg, A, B, C);

    sk=0;
    printf("\nC=\n");
    for(i=0;i<Deg;i++){
        for(j=0;j<Deg;j++){
            printf("%lg,", C[sk++]);
        }
        printf("\n");
    }

    free(A);
    free(B);
    free(C);
}

void inverseMatrix(void){
    doublereal *Matrix;
    int i,j, sk=0;
    Matrix = (doublereal*) malloc(sizeof(doublereal)*Deg*Deg);

    for(i=0;i<Deg;i++){
        for(j=0;j<Deg;j++){
            Matrix[sk] = sk;
            sk++;
        }
        Matrix[sk-1] = sk;
    }
    sk = 0;
    for(i=0;i<Deg;i++){
        for(j=0;j<Deg;j++){
            Matrix[sk] *= Matrix[sk];
            printf("%f, ", Matrix[sk++]);
        }
        printf("\n");
    }
    invM(Matrix, Deg);


//    printf("LWORK[%d]=%f")
    sk=0;
    for(i=0;i<Deg;i++){
        for(j=0;j<Deg;j++){
            printf("%lg,", Matrix[sk++]);
        }
        printf("\n");
    }


//    free(Matrix);
    free(Matrix);

}

#define Np (Deg+1)*(Deg+2)/2

void GetTriCoordTest(void){
//    int N = 5;
    double *r = BuildVector(Np);
    double *s = BuildVector(Np);
    int i;


    GetTriCoord(Deg, r, s);

    for(i=0;i<Np;i++){
        printf("r[%d]=%12.4f, s[%d]=%12.4f\n", i,r[i],i,s[i]);
    }

    DestroyVector(r);
    DestroyVector(s);

}

void WarpfactorTest(void){
    int i;
    double *v = BuildVector(Np);
    double *w = BuildVector(Np);

    for(i=0;i<Np;i++){
        v[i] = i/(Np-1.0)*2.0 - 1.0;
        printf("v[%d] = %lg,\t w[%d] = %lg\n", i, v[i], i, w[i]);
    }

    Warpfactor(Deg, v, Np, w);

    for(i=0;i<Np;i++){
        printf("v[%d] = %lg,\t w[%d] = %e\n", i, v[i], i, w[i]);
    }

    DestroyVector(v);
    DestroyVector(w);
}
