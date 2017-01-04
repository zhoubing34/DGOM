//
// Created by li12242 on 12/19/16.
//

#include <StandCell/sc_stdcell.h>
#include "mr_mesh.h"
#include "mr_mesh_cellConnect.h"
#include "mr_mesh_nodeConnect.h"

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

    /* initialize boundary condition */
    mesh->EToBS = NULL;
    mesh->bcIndList = NULL;
    mesh->obvertlist = NULL;

    return mesh;
}

void mr_mesh_free(parallMesh *mesh){
    IntMatrix_free(mesh->EToE);
    IntMatrix_free(mesh->EToF);
    IntMatrix_free(mesh->EToP);

    IntVector_free(mesh->Npar);
    IntVector_free(mesh->vmapM);
    IntVector_free(mesh->vmapP);

    IntVector_free(mesh->nodeIndexOut);
    IntVector_free(mesh->cellIndexIn);
    IntVector_free(mesh->cellIndexOut);

    return;
}