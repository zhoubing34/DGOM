#ifndef UNSTRCUTMESH_H
#define UNSTRCUTMESH_H

#include "StandCell/sc_stdcell.h"
#include "mpi.h"

/* geometry grid structure */
typedef struct{
    int dim; ///> dimensions
    cellType type; ///> type of element
    int nprocs; ///> number of process
    int procid; ///> index of process
    stdCell *cell; ///> standard element
    int Nv;///> Num of vertex
    int K; ///> Num of element
    int **EToV; ///> vertex list in elements, the index start from 0
    double *vx, *vy, *vz; ///> coordinate of vertex
} geoGrid;


geoGrid* mr_grid_create(stdCell *shape, int K, int Nv, double *vx, double *vy, double *vz, int **EToV);
void mr_grid_free(geoGrid *grid);

#endif // UNSTRCUTMESH_H
