#ifndef DGOM_PHYSDOMAIN_H
#define DGOM_PHYSDOMAIN_H

#include "MultiRegions/Edge/dg_edge.h"


typedef struct {
    dg_real *px_Q; ///< dfdx partial derivative for x direction
    dg_real *py_Q; ///< dfdy partial derivative for y direction
    dg_real *px_inQ, *px_outQ; ///< send and recv buffers for p_Q
    dg_real *py_inQ, *py_outQ; ///< send and recv buffers for q_Q
    dg_real *vis_Q; ///< viscosity on each node
} pf_LDG_solver;

/**
 * @brief physical field structure.
 */
typedef struct dg_phys{
    int Nfield; ///< number of variable fields;

    dg_edge *edge; ///< parallel edge structure;
    dg_mesh *mesh; ///< parallel mesh structure;
    dg_region *region; ///< multi-region structure;
    dg_grid *grid; ///< geometry grid structure;
    dg_cell *cell; ///< standard cell structure;

    /* cell info */
    dg_real *c_Q; ///< cell information;
    dg_real *c_recvQ; ///< recv buffers for cell information;

    dg_real *f_Q; ///< nodal information;
    dg_real *f_recvQ; ///< recv buffers for nodal information;
    dg_real *f_extQ; ///< external nodal data;

    dg_real *f_rhsQ; ///< RHS data;
    dg_real *f_resQ; ///< residual data;

    /** function to fetch node buffer with other process; */
    int (*fetch_node_buffer)(struct dg_phys *phys, MPI_Request *send_requests, MPI_Request *recv_requests);
    /** function to fetch cell buffer with other process; */
    int (*fetch_cell_buffer)(struct dg_phys *phys, MPI_Request *send_requests, MPI_Request *recv_requests);

} dg_phys;

dg_phys* dg_phys_create(int Nfields, dg_edge *edge);
void dg_phys_free(dg_phys *phys);

#define dg_phys_Nfield(phys) phys->Nfield
#define dg_phys_f_Q(phys) phys->f_Q
#define dg_phys_f_extQ(phys) phys->f_extQ
#define dg_phys_f_recvQ(phys) phys->f_recvQ
#define dg_phys_f_rhsQ(phys) phys->f_rhsQ
#define dg_phys_c_Q(phys) phys->c_Q
#define dg_phys_c_recvQ(phys) phys->c_recvQ

#endif //DGOM_PHYSDOMAIN_H
