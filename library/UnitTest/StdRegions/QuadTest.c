//
// Created by li12242 on 16/6/30.
//

#include "QuadTest.h"

#define RETURN printf("\n");

int main(void){
    StdRegions2d *quad = StdQuadEle_create(Deg);

    /* parameters */

    /* Fmask */
    int Nfaces  = Dnfaces, Nfp = Dfp;
    int **Fmask = IntMatrix_create(Nfaces, Nfp);

    int i,j;
    for(i=0;i<Nfaces;i++){
        for(j=0;j<Nfp;j++)
            Fmask[i][j] = TestQuad_fmask[i][j];
    }
    CreateIntMatrixTest("Quad Fmask", quad->Fmask, Fmask, Nfaces, Nfp); RETURN
    IntMatrix_free(Fmask);

    /*Coordinate*/
    int    Np = (Deg+1)*(Deg+1);
    double *r = Vector_create(Np);
    double *s = Vector_create(Np);

    for(i=0;i<Np;i++) {
        r[i] = TestQuad_r[i];
        s[i] = TestQuad_s[i];
    }
    CreateVectorTest("Coordinate r", quad->r, r, Dnp); RETURN
    CreateVectorTest("Coordinate s", quad->s, s, Dnp); RETURN
    Vector_free(r);
    Vector_free(s);

    /* matrix */
    double **V  = Matrix_create(Np, Np);
    double **M  = Matrix_create(Np, Np);
    double **Dr = Matrix_create(Np, Np);
    double **Ds = Matrix_create(Np, Np);

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

    Matrix_free(V);
    Matrix_free(M);
    Matrix_free(Dr);
    Matrix_free(Ds);

    /* LIFT test */
    double **LIFT = Matrix_create(Np, Nfp * Nfaces);

    for(i=0;i<Np;i++){
        for(j=0;j<Nfp*Nfaces;j++){
            LIFT[i][j] = TestQuad_LIFT[i][j];
        }
    }
    CreateMatrixTest("LIFT matrix", quad->LIFT, LIFT, Dnp, Dfp*Dnfaces); RETURN
    Matrix_free(LIFT);

    /* finish test */
    StdRegions2d_free(quad);

    return 0;
}