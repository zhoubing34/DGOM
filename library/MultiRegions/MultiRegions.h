/**
 * @file
 * Mesh object for multi-process
 *
 * @brief
 * The basic mesh information
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#ifndef DGOM_MULTIREGIONS_H
#define DGOM_MULTIREGIONS_H

#include "StdRegions/StdRegions.h"
#include "LocalRegions/LocalRegions.h"
#include <mpi.h>
#include <parmetisbin.h>

typedef struct MultiReg2d{

    /** index of local process */
    int procid;
    /** number of processes */
    int nprocs;

    /** standard element */
    StdRegions2d *StdElement;
    /** total number of vertex */
    int Nv;
    /** total number of elements */
    int K;

    /* coordinate */

    /** vertex coordinate */
    double **GX;
    /** vertex coordinate */
    double **GY;
    /** node coordinate */
    double **x;
    /** node coordinate */
    double **y;

    /* mesh and nodes info */

    /* element connection info */
    int **EToE, **EToV, **EToP;

    /** number of faces to send recv to each proc */
    int *Npar;
    /** local element index of parallel nodes */
    int **parK;
    /** local face index of parallel nodes */
    int **parF;

    /** total number of nodes to send recv */
    int parNtotalout;
    /** index list of nodes to send out */
    int *parmapOUT;

    /** volume id of -ve trace of face node */
    int *vmapM;
    /** volume id of +ve trace of face node */
    int *vmapP;

    /* info of each element */

    /** jacobi on each nodes */
    double *J;
    /** radius of the circumscribed circle */
    float *ciradius;
    /* area of each ele */
    float *area;


}MultiReg2d;

/* public function */
MultiReg2d* GenMultiReg2d
        (StdRegions2d *shape, int K, int Nv, int **EToV, double *VX, double *VY);

void FreeMultiReg2d(MultiReg2d *mesh);
void LoadBalance(MultiReg2d *mesh);

#endif //DGOM_MULTIREGIONS_H
