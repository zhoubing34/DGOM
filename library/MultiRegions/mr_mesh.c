//
// Created by li12242 on 12/19/16.
//

#include <StandCell/sc_stdcell.h>
#include "mr_mesh.h"
#include "mr_mesh_cellConnect.h"
#include "mr_mesh_nodeConnect.h"
#include "mr_grid.h"

/**
 * @brief create mesh object
 * @param[in] region multi-region object
 */
parallMesh* mr_mesh_create(multiReg *region){

    parallMesh *mesh = calloc(1, sizeof(parallMesh));

    mesh->dim = region->dim;
    mesh->nprocs = region->nprocs;
    mesh->procid = region->procid;

    /* objects */
    mesh->region = region;
    mesh->grid = region->grid;
    mesh->cell = region->cell;

    /* element pairs, including the adjacent element in other process */
    mr_mesh_cellConnect2d(mesh);

    /* node pairs, including the adjacent nodes in other process */
    mr_mesh_nodeConnect2d(mesh);

    return mesh;
}