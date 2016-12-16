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

#include "StandCell/sc_stdcell.h"
#include "LocalRegions/LocalRegions.h"
#include <mpi.h>
#include "UnstructMesh.h"

#define NODETOL  1e-6f

typedef struct {

    /** index of local process */
    int procid;
    /** number of processes */
    int nprocs;

    /** standard element */
    stdCell *stdcell;
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

    /* info of each element */

    /** jacobi on each nodes */
    double *J;
    /** radius of the circumscribed circle */
    double *ciradius;
    /* area of each ele */
    double *area;
    /* volume geometric factors */
    int Nvgeo;
    real *vgeo;

}MultiReg2d;

/* public function */
MultiReg2d* MultiReg2d_create(stdCell *shape, UnstructMesh *grid);
void MultiReg2d_free(MultiReg2d *mesh);

#endif //DGOM_MULTIREGIONS_H
