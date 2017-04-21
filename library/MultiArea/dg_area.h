//
// Created by li12242 on 17/4/17.
//

#ifndef DGOM_DG_AREA_H
#define DGOM_DG_AREA_H

#include "Grid/dg_grid.h"
#include "Region/dg_region.h"
#include "Mesh/dg_mesh.h"
#include "Edge/dg_edge.h"

typedef struct dg_area{
    int procid; ///< process id;
    int nprocs; ///< number of processes;

    dg_cell *cell;
    dg_grid *grid;
    dg_region *region;
    dg_mesh *mesh;
    dg_edge *edge;
} dg_area;

dg_area *dg_area_create_from_file(dg_cell *cell, char *casename);
dg_area *dg_area_create_uniform(dg_cell *cell, int Mx, int My,
                                double xmin, double xmax,
                                double ymin, double ymax, int type);
void dg_area_free(dg_area *area);

#define dg_area_grid(area) area->grid
#define dg_area_cell(area) area->cell
#define dg_area_region(area) area->region
#define dg_area_mesh(area) area->mesh
#define dg_area_edge(area) area->edge
#define dg_area_nprocs(area) area->nprocs
#define dg_area_procid(area) area->procid

#endif //DGOM_DG_AREA_H
