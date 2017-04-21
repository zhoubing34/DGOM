//
// Created by li12242 on 17/4/17.
//

#include "dg_area.h"

/**
 * @brief create dg_area strcture from input file
 * @param[in] cell pointer to dg_cell structure;
 * @param[in] casename name of input case;
 * @return area pointer to a new dg_area structure.
 */
dg_area * dg_area_create_from_file(dg_cell *cell, char *casename){
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    dg_area *area = (dg_area*) calloc(1, sizeof(dg_area));

    area->cell = cell;
    area->procid = procid;
    area->nprocs = nprocs;

    dg_grid *grid = NULL;
    switch (dg_cell_celltype(cell)){
        case TRIANGLE:
            grid = dg_grid_create_from_file2d(cell, casename); break;
        case QUADRIL:
            grid = dg_grid_create_from_file2d(cell, casename); break;
        default:
            fprintf(stderr, "%s (%d): Unknown cell type %d\n",
                    __FUNCTION__, __LINE__, dg_cell_celltype(cell));
            exit(-1);
    }
    dg_region *region = dg_region_create(grid);
    dg_mesh *mesh = dg_mesh_create(region);
    dg_edge *edge = dg_edge_create(mesh);

    area->grid = grid;
    area->region = region;
    area->mesh = mesh;
    area->edge = edge;
    return area;
}

dg_area * dg_area_create_uniform(dg_cell *cell, int Mx, int My,
                                 double xmin, double xmax,
                                 double ymin, double ymax,
                                 int type){
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    dg_area *area = (dg_area*) calloc(1, sizeof(dg_area));

    area->cell = cell;
    area->procid = procid;
    area->nprocs = nprocs;

    dg_grid *grid;
    switch (dg_cell_celltype(cell)){
        case TRIANGLE:
            grid = dg_grid_create_uniform_tri(cell, Mx, My, xmin, xmax, ymin, ymax, type); break;
        case QUADRIL:
            grid = dg_grid_create_uniform_quad(cell, Mx, My, xmin, xmax, ymin, ymax); break;
        default:
            fprintf(stderr, "%s (%d): Unknown cell type %d\n",
                    __FUNCTION__, __LINE__, dg_cell_celltype(cell));
            exit(-1);
    }

    dg_region *region = dg_region_create(grid);
    dg_mesh *mesh = dg_mesh_create(region);
    dg_edge *edge = dg_edge_create(mesh);

    area->grid = grid;
    area->region = region;
    area->mesh = mesh;
    area->edge = edge;
    return area;
}


void dg_area_free(dg_area *area){
    dg_grid_free(area->grid);
    dg_region_free(area->region);
    dg_mesh_free(area->mesh);
    dg_edge_free(area->edge);
    free(area);
    return;
}