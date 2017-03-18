/**
 * @file
 * physical fields
 * @brief
 *
 * @author li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "dg_phys.h"
#include "dg_phys_fetch_buffer.h"

#define DEBUG 0

dg_phys* dg_phys_create(int Nfields, dg_edge *edge){

    dg_phys *phys = (dg_phys *) calloc(1, sizeof(dg_phys));

    phys->edge = edge;
    phys->mesh = edge->mesh;
    phys->region = edge->region;
    phys->grid = edge->grid;
    phys->cell = edge->cell;

    phys->Nfield = Nfields; ///< number of fields

    const int K = dg_grid_K(phys->grid);
    const int Np = dg_cell_Np(phys->cell);
    const int Nparface = dg_mesh_NfetchFace(phys->mesh);
    const int Nparnode = dg_mesh_NfetchNode(phys->mesh);

    /* nodal array */
    phys->f_Q    = (dg_real *) calloc((size_t) K*Np*Nfields, sizeof(dg_real));
    phys->f_rhsQ = (dg_real *) calloc((size_t) K*Np*Nfields, sizeof(dg_real));
    phys->f_resQ = (dg_real *) calloc((size_t) K*Np*Nfields, sizeof(dg_real));
    phys->f_extQ = (dg_real *) calloc((size_t) K*Np*Nfields, sizeof(dg_real));
    phys->f_recvQ = (dg_real *) calloc((size_t) Nparnode*Nfields, sizeof(dg_real));

    /* volume array */
    phys->c_Q = (dg_real *) calloc((size_t) K*Nfields, sizeof(dg_real));
    phys->c_recvQ = (dg_real *) calloc((size_t) Nparface*Nfields, sizeof(dg_real));

    /* fetch buffer */
    phys->fetch_node_buffer = dg_phys_fetch_node_buffer;
    phys->fetch_cell_buffer = dg_phys_fetch_cell_buffer;
    return phys;
}

void dg_phys_free(dg_phys *phys){

    free(phys->f_Q);
    free(phys->f_rhsQ);
    free(phys->f_resQ);
    free(phys->f_extQ);
    free(phys->f_recvQ);

    free(phys->c_Q);
    free(phys->c_recvQ);
    free(phys);
    return;
}
