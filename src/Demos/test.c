
#include <stdio.h>
#include "LibUtilities/LibUtilities.h"
#include "StdRegions/StdRegions.h"

#define Deg 2

int main(int argc, char **argv) {

    void TriEleTest(void);
    TriEleTest();

    return 0;
}


void TriEleTest(void){

    StdRegions2d *tri;
    tri = GenStdTriEle(Deg);

    PrintMatrix("Mass matrix", tri->M, tri->Np, tri->Np);
    printf("\n\n");
    PrintMatrix("Dr", tri->Dr, tri->Np, tri->Np);
    printf("\n\n");
    PrintMatrix("Ds", tri->Ds, tri->Np, tri->Np);
    printf("\n\n");

    int i,j;
    printf("Fmask = \n");
    for(i=0;i<tri->Nfaces;i++){
        for(j=0;j<tri->Nfp;j++){
            printf("%d,", tri->Fmask[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");

    printf("\n\n");
    PrintMatrix("LIFT", tri->LIFT, tri->Np, tri->Nfp*tri->Nfaces);

    FreeStdRegions2d(tri);
}

void MatrixMuti(void){
    int i,j,sk=0;
    const int M=Deg;
    double *A, *B, *C;
    A = (double*) malloc(sizeof(double)*Deg*Deg);
    B = (double*) malloc(sizeof(double)*Deg*Deg);
    C = (double*) malloc(sizeof(double)*Deg*Deg);

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
    double *Matrix;
    int i,j, sk=0;
    Matrix = (double*) malloc(sizeof(double)*Deg*Deg);

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

//void GetTriCoordTest(void){
////    int N = 5;
//    double *r = BuildVector(Np);
//    double *s = BuildVector(Np);
//    int i;
//
//
//    GetTriCoord(Deg, r, s);
//
//    for(i=0;i<Np;i++){
//        printf("r[%d]=%12.4f, s[%d]=%12.4f\n", i,r[i],i,s[i]);
//    }
//
//    DestroyVector(r);
//    DestroyVector(s);
//
//}
//
//void WarpfactorTest(void){
//    int i;
//    double *v = BuildVector(Np);
//    double *w = BuildVector(Np);
//
//    for(i=0;i<Np;i++){
//        v[i] = i/(Np-1.0)*2.0 - 1.0;
//        printf("v[%d] = %lg,\t w[%d] = %lg\n", i, v[i], i, w[i]);
//    }
//
//    Warpfactor(Deg, v, Np, w);
//
//    for(i=0;i<Np;i++){
//        printf("v[%d] = %lg,\t w[%d] = %e\n", i, v[i], i, w[i]);
//    }
//
//    DestroyVector(v);
//    DestroyVector(w);
//}
