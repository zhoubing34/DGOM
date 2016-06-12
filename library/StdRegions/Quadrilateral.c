/**
 * @file
 * Quadrilateral.c
 *
 * @brief
 * Functions related with standard quadrilateral elements
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "StdRegions.h"

/**
 * @brief
 * Generation of standard triangle element
 *
 * @details
 *
 * @param[in] N polynomial order
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * tri | StdRegions2d* |
 *
 */
StdRegions2d* GenStdQuadEle(int N){
    StdRegions2d *quad = (StdRegions2d *) calloc(1, sizeof(StdRegions2d));

    int Np, Nfaces, Nv, Nfp;
    /* basic info */
    quad->N = N;
    Np = (N+1)*(N+1);
    quad->Nfaces = 4;
    quad->Nfp = N+1;
    quad->Nv = 4;

    /* coordinate */
    quad->r = (double *) calloc(Np, sizeof(double));
    quad->s = (double *) calloc(Np, sizeof(double));
    return quad;
}