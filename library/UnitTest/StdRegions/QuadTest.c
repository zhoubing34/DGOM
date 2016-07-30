//
// Created by li12242 on 16/6/30.
//

#include "QuadTest.h"

#define RETURN printf("\n");

int main(void){
    StdRegions2d *quad = GenStdQuadEle(Deg);

    /* parameters */

    /* Fmask */
    int Nfaces  = Dnfaces, Nfp = Dfp;
    int **Fmask = BuildIntMatrix(Nfaces, Nfp);

    int i,j;
    for(i=0;i<Nfaces;i++){
        for(j=0;j<Nfp;j++)
            Fmask[i][j] = TestQuad_fmask[i][j];
    }
    CreateIntMatrixTest("Quad Fmask", quad->Fmask, Fmask, Nfaces, Nfp); RETURN
    DestroyIntMatrix(Fmask);

    /*Coordinate*/
    int    Np = (Deg+1)*(Deg+1);
    double *r = BuildVector(Np);
    double *s = BuildVector(Np);

    for(i=0;i<Np;i++) {
        r[i] = TestQuad_r[i];
        s[i] = TestQuad_s[i];
    }
    CreateVectorTest("Coordinate r", quad->r, r, Dnp); RETURN
    CreateVectorTest("Coordinate s", quad->s, s, Dnp); RETURN
    DestroyVector(r);
    DestroyVector(s);

    /* matrix */
    double **V  = BuildMatrix(Np, Np);
    double **M  = BuildMatrix(Np, Np);
    double **Dr = BuildMatrix(Np, Np);
    double **Ds = BuildMatrix(Np, Np);

    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            V [i][j] = TestQuad_V[i][j];
            M [i][j] = TestQuad_M[i][j];
            Dr[i][j] = TestQuad_Dr[i][j];
            Ds[i][j] = TestQuad_Ds[i][j];
        }
    }
    CreateMatrixTest("Vandermonde matrix", quad->V, V, Dnp, Dnp); RETURN
    CreateMatrixTest("Vandermonde matrix", quad->M, M, Dnp, Dnp); RETURN
    CreateMatrixTest("Drivative matrix Dr", quad->Dr, Dr, Dnp, Dnp); RETURN
    CreateMatrixTest("Drivative matrix Ds", quad->Ds, Ds, Dnp, Dnp); RETURN

    DestroyMatrix(V);
    DestroyMatrix(M);
    DestroyMatrix(Dr);
    DestroyMatrix(Ds);

    /* LIFT test */
    double **LIFT = BuildMatrix(Np, Nfp*Nfaces);

    for(i=0;i<Np;i++){
        for(j=0;j<Nfp*Nfaces;j++){
            LIFT[i][j] = TestQuad_LIFT[i][j];
        }
    }
    CreateMatrixTest("LIFT matrix", quad->LIFT, LIFT, Dnp, Dfp*Dnfaces); RETURN
    DestroyMatrix(LIFT);

    /* finish test */
    FreeStdRegions2d(quad);

    return 0;
}