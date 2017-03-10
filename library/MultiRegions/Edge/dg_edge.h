//
// Created by li12242 on 17/3/10.
//

#ifndef DGOM_DG_EDGE_H
#define DGOM_DG_EDGE_H

#include "MultiRegions/Mesh/dg_mesh.h"

typedef struct dg_edge{
    int Nedge; ///< total number of edges
    int Nnode; ///< total number of nodes
    int procid;
    int nprocs;

    dg_cell *cell;
    dg_grid *grid;
    dg_region *region;
    dg_mesh *mesh;

    int *varkM, *varkP; ///< cell index of local and adjacent cell
    int *varfM, *varfP; ///< face index of local and adjacent cell
    int *varpM, *varpP; ///< node index of local and adjacent node
    int *ftype; ///< face type of each edges
    dg_real *fsc; ///< face Jacobi divided by volume Jacobi
    dg_real *nx, *ny, *nz; ///< outward normal vector of local cell

    void(*free_func)(struct dg_edge *edge);
}dg_edge;

dg_edge *dg_edge_create(dg_mesh *mesh);
void dg_edge_free(dg_edge *edge);

#define dg_edge_Nedge(edge) edge->Nedge
#define dg_edge_Nnode(edge) edge->Nnode
#define dg_edge_procid(edge) edge->procid
#define dg_edge_nprocs(edge) edge->nprocs


#endif //DGOM_DG_EDGE_H
