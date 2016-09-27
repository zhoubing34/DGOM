#include "StdRegions.h"

/**
 * @brief
 * Allocate stdRegions
 */
StdRegions2d* StdRegions2d_create(int N, ElementType eleType){
    StdRegions2d *std;
    switch (eleType){
        case TRIANGLE:
            std = StdTriEle_create(N); break;
        case QUADRIL:
            std = StdQuadEle_create(N); break;
        default:
            printf("Unknown element type\n");
            exit(-1);
    }
    return std;
}

/**
 * @brief
 * Free variables fields in StdRegions2d
 */
void StdRegions2d_free(StdRegions2d *triangle){
    /* fmaks */
    IntMatrix_free(triangle->Fmask);
    /* coordinate */
    Vector_free(triangle->r);
    Vector_free(triangle->s);
    /* vandermonde matrix */
    Matrix_free(triangle->V);
    /* mass matrix */
    Matrix_free(triangle->M);
    /* Derivative Matrix */
    Matrix_free(triangle->Dr);
    Matrix_free(triangle->Ds);

    /* LIFT */
    Matrix_free(triangle->LIFT);

    /* float version */
    free(triangle->f_LIFT);
    free(triangle->f_Dr);
    free(triangle->f_Ds);

}