#ifndef UNSTRCUTMESH_H
#define UNSTRCUTMESH_H

#include "StandCell/dg_cell.h"
#include "mpi.h"

/* geometry grid structure */
typedef struct{
    int dim; ///> dimensions
    //dg_cell_type type; ///> type of element
    int nprocs; ///> number of process
    int procid; ///> index of process
    dg_cell *cell; ///> standard element
    int Nv;///> Num of vertex
    int K; ///> Num of element
    int **EToV; ///> vertex list in elements, the index start from 0
    double *vx, *vy, *vz; ///> coordinate of vertex
} dg_grid;


dg_grid* mr_grid_create(dg_cell *shape, int K, int Nv, double *vx, double *vy, double *vz, int **EToV);
void mr_grid_free(dg_grid *grid);

#endif // UNSTRCUTMESH_H
