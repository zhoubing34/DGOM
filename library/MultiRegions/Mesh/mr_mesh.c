//
// Created by li12242 on 12/19/16.
//

#include <StandCell/dg_cell.h>
#include "mr_mesh.h"
#include "mr_mesh_cellConnect.h"
#include "mr_mesh_nodeConnect.h"

/**
 * @brief create mesh object
 * @param[in] region multi-region object
 */
dg_mesh* mr_mesh_create(dg_region *region){

    dg_mesh *mesh = calloc(1, sizeof(dg_mesh));

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

    /* initialize boundary condition */
    mesh->EToBS = NULL;
    mesh->bcind = NULL;
    mesh->obcind = NULL;

    return mesh;
}

void mr_mesh_free(dg_mesh *mesh){
    matrix_int_free(mesh->EToE);
    matrix_int_free(mesh->EToF);
    matrix_int_free(mesh->EToP);

    vector_int_free(mesh->Npar);
    vector_int_free(mesh->vmapM);
    vector_int_free(mesh->vmapP);

    vector_int_free(mesh->nodeIndexOut);
    vector_int_free(mesh->cellIndexIn);
    vector_int_free(mesh->cellIndexOut);

    return;
}