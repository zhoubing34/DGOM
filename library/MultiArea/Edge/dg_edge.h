//
// Created by li12242 on 17/3/10.
//

#ifndef DGOM_DG_EDGE_H
#define DGOM_DG_EDGE_H

#include "../Mesh/dg_mesh.h"

/**
 * @brief
 */
typedef struct dg_edge{
    int Nedge; ///< total number of edges;
    int Nnode; ///< total number of nodes;
//    int procid; ///< process id;
//    int nprocs; ///< number of process;

    dg_mesh *mesh;
    /* face info */
    int *varkM,*varkP; ///< cell index of local and adjacent cell;
    int *varfM,*varfP; ///< face index of local and adjacent cell;
    int *ftype; ///< face type of each edges;
    int *surfinfo; ///< gather the varkM/varkP, varfM/varfP and ftype;
    /* node info */
    int *varpM,*varpP; ///< node index of local and adjacent node;
    int *varfpM,*varfpP; ///< face node index of local and adjacent node;
    dg_real *fsc; ///< face Jacobi divided by volume Jacobi;
    dg_real *nx,*ny,*nz; ///< outward normal vector of local cell;
    dg_real *nodeinfo; ///< gather the varpM/varpP, fsc and nx/ny/nz;

}dg_edge;

dg_edge *dg_edge_create(dg_mesh *mesh);
void dg_edge_free(dg_edge *edge);

#define dg_edge_mesh(edge) edge->mesh
#define dg_edge_region(edge) dg_mesh_region(dg_edge_mesh(edge))
#define dg_edge_grid(edge) dg_mesh_grid(dg_edge_mesh(edge))
#define dg_edge_cell(edge) dg_mesh_cell(dg_edge_mesh(edge))
#define dg_edge_Nedge(edge) edge->Nedge
#define dg_edge_Nnode(edge) edge->Nnode
#define dg_edge_procid(edge) edge->procid
#define dg_edge_nprocs(edge) edge->nprocs


#endif //DGOM_DG_EDGE_H
