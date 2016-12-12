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
#include "UnstructMesh.h"

#define NODETOL  1e-6f

typedef struct {

    /** index of local process */
    int procid;
    /** number of processes */
    int nprocs;

    /** standard element */
    StdRegions2d *stdcell;
    /** total number of vertex */
    int Nv;
    /** local number of elements */
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
    int **EToV;
    int **EToE;
    int **EToF;
    int **EToP;

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

    /* total number of element's value to send and recv */
    int parEtotalout;
    /* index list of elements to send out */
    int *elemapOut;

    /** volume id of -ve trace of face node */
    int *vmapM;
    /** volume id of +ve trace of face node */
    int *vmapP;

    /* info of each element */

    /** jacobi on each nodes */
    double *J;
    /** radius of the circumscribed circle */
    double *ciradius;
    /* area of each ele */
    double *area;
    /* volume geometric factors */
    real *vgeo;

}MultiReg2d;



/* public function */
MultiReg2d* MultiReg2d_create(StdRegions2d *shape, UnstructMesh *grid);
void MultiReg2d_free(MultiReg2d *mesh);

#endif //DGOM_MULTIREGIONS_H
