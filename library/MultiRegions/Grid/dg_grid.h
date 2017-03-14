#ifndef UNSTRCUTMESH_H
#define UNSTRCUTMESH_H

#include "StandCell/dg_cell.h"
#include "mpi.h"

typedef enum {
    INNERLOC=0, // 0-inner face2d, adjacent face2d in the same process
    INNERBS=1,  // 1-inner boundary surface, adjacent face2d in other process
    SLIPWALL=2, // 2-slip wall
    NSLIPWALL=3, // 3-non-slip wall
    OPENBS = 4   // 4-open boundary
} dg_face_type;

typedef enum{
    COMPUTE_CELL = 0,
    SPONGE_CELL = 1,
    REFINED_CELL = 2,
    TROUBLE_CELL = 3
} dg_grid_cell_type;

/* geometry grid structure */
typedef struct dg_grid{
    dg_cell *cell; ///< standard element

    int nprocs; ///< number of process
    int procid; ///< index of process
    int Nv;///< Num of vertex
    int K; ///< Num of element
    int **EToV; ///< vertex list in elements, start from 0
    int **EToE; ///< adjacent element id
    int **EToF; ///< adjacent face id (0,1,2...)
    int **EToP; ///< process id of adjacent element
    int **EToBS; ///< boundary surface type of adjacent element
    int *EToR; ///< region id of each cell
    double *vx, *vy, *vz; ///< vertex coordinate

    void (*add_BS)(struct dg_grid *grid, int Nsurf, int **SFToV);
    void (*free_func)(struct dg_grid *grid); ///< free function
} dg_grid;


dg_grid* dg_grid_create(dg_cell *cell, int K, int Nv, double *vx, double *vy, double *vz, int **EToV);
void dg_grid_add_BS(dg_grid *grid, int Nsurf, int **SFToBS);
void dg_grid_free(dg_grid *grid);

#define dg_grid_K(grid) grid->K
#define dg_grid_Nv(grid) grid->Nv
#define dg_grid_procid(grid) grid->procid
#define dg_grid_nprocs(grid) grid->nprocs

#endif // UNSTRCUTMESH_H
