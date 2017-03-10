//
// Created by li12242 on 12/19/16.
//

#ifndef DGOM_MR_MESH_H
#define DGOM_MR_MESH_H

#include "MultiRegions/Region/dg_region.h"

typedef struct dg_mesh{

    int procid; ///< process id
    int nprocs; ///< number of process

    dg_region *region; ///< multi-region object
    dg_grid *grid; ///< geometry grid object
    dg_cell *cell; ///< standard element object

    int *parfaceNum; ///< number of faces adjacent to each process
    int TotalParFace; ///< total number of parallel faces
    int *parcell; ///< map of cell index for buffer
    int *parface; ///< map of local face (0,1,2..) index for buffer

    int *parnodeNum; ///< number of nodes adjacent to each process
    int TotalParNode; ///< total number of nodes to send recv
    int *parnode; ///< index list of nodes to send out

    void (*free_func)(struct dg_mesh *mesh);
} dg_mesh;

dg_mesh* dg_mesh_create(dg_region *region);
void dg_mesh_free(dg_mesh *mesh);

#define dg_mesh_procid(mesh) mesh->procid
#define dg_mesh_nprocs(mesh) mesh->nprocs
#define dg_mesh_Nparf(mesh) mesh->TotalParFace
#define dg_mesh_Nparn(mesh) mesh->TotalParNode

#endif //DGOM_MR_MESH_H
