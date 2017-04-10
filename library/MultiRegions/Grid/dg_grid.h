#ifndef UNSTRCUTMESH_H
#define UNSTRCUTMESH_H

#include "StandCell/dg_cell.h"
#include "mpi.h"

/**
 * @brief cell surface type list.
 */
typedef enum {
    FACE_INNER = 0, ///< inner local faces, adjacent to other cells;
    FACE_PARALL = 1, ///< inner boundary surface, adjacent to other processes;
    FACE_SLIPWALL = 2, ///< slip wall;
    FACE_NSLIPWALL = 3, ///< non-slip wall;
    FACE_OPENBS = 4, ///< open boundary;
} dg_grid_face_type;

/**
 * @brief cell type list.
 */
typedef enum{
    CELL_COMPUTE = 0, ///< normal compute cell;
    CELL_SPONGE = 1, ///< sponge region cell;
    CELL_REFINED = 2, ///< refined cell;
    CELL_TROUBLE = 3, ///< tourble cell;
} dg_grid_cell_type;

/**
 * @brief geometry grid structure
 */
typedef struct dg_grid{
    dg_cell *cell; ///< standard element;

    int nprocs; ///< number of process;
    int procid; ///< index of process;
    int Nv;///< Num of vertex;
    int K; ///< Num of element;
    int **EToV; ///< vertex list in elements, start from 0;
    int **EToE; ///< adjacent element id;
    int **EToF; ///< adjacent face id (0,1,2...);
    int **EToP; ///< process id of adjacent element;
    int **EToBS; ///< boundary surface type of adjacent element;
    int *EToR; ///< region id of each cell;
    double *vx, *vy, *vz; ///< vertex coordinate;

    /** add boundary surface condition */
    void (*add_BS)(struct dg_grid *grid, int Nsurface, int **SFToV);
    /** add boundary surface condition from input file */
    void (*add_BS_from_file)(struct dg_grid *grid, char *casename);
    /** project the grid vertex values to nodes */
    void (*proj_vert2node)(struct dg_grid *grid, int Nfield, double *vertval, double *nodeval);
} dg_grid;


dg_grid* dg_grid_create(dg_cell *cell, int K, int Nv, double *vx, double *vy, double *vz, int **EToV);
void dg_grid_free(dg_grid *grid);

#include "dg_grid_reader.h"

#define dg_grid_cell(grid) grid->cell
#define dg_grid_K(grid) grid->K
#define dg_grid_Nv(grid) grid->Nv
#define dg_grid_EToV(grid) grid->EToV
#define dg_grid_EToE(grid) grid->EToE
#define dg_grid_EToF(grid) grid->EToF
#define dg_grid_EToP(grid) grid->EToP
#define dg_grid_EToBS(grid) grid->EToBS
#define dg_grid_EToR(grid) grid->EToR

#define dg_grid_procid(grid) grid->procid
#define dg_grid_nprocs(grid) grid->nprocs
#define dg_grid_vx(grid) grid->vx
#define dg_grid_vy(grid) grid->vy
#define dg_grid_vz(grid) grid->vz

#endif // UNSTRCUTMESH_H
