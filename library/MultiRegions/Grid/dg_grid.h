#ifndef UNSTRCUTMESH_H
#define UNSTRCUTMESH_H

#include "StandCell/dg_cell.h"
#include "mpi.h"

/* geometry grid structure */
typedef struct dg_grid{
    dg_cell *cell; ///< standard element

    int nprocs; ///< number of process
    int procid; ///< index of process
    int Nv;///< Num of vertex
    int K; ///< Num of element
    int **EToV; ///< vertex list in elements, start from 0
    double *vx, *vy, *vz; ///< vertex coordinate

    void (*free_func)(struct dg_grid *grid); ///< free function
} dg_grid;


dg_grid* dg_grid_create(dg_cell *cell, int K, int Nv, double *vx, double *vy, double *vz, int **EToV);
void dg_grid_free(dg_grid *grid);

#endif // UNSTRCUTMESH_H
