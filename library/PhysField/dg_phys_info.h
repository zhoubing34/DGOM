//
// Created by li12242 on 17/3/22.
//

#ifndef DGOM_DG_PHYS_INFO_H
#define DGOM_DG_PHYS_INFO_H

#include "MultiArea/dg_area.h"

typedef struct dg_phys_info{
    int Nfield;

//    dg_edge *edge; ///< parallel edge structure;
//    dg_mesh *mesh; ///< parallel mesh structure;
//    dg_region *region; ///< multi-region structure;
//    dg_grid *grid; ///< geometry grid structure;
//    dg_cell *cell; ///< standard cell structure;
    dg_area *area;

    dg_real *c_Q; ///< cell information;
    dg_real *c_recvQ; ///< recv buffers for cell information;

    dg_real *f_Q; ///< nodal information;
    dg_real *f_recvQ; ///< recv buffers for nodal information;

    dg_real *f_rhsQ; ///< right hand side data;
    dg_real *f_resQ; ///< residual data;

    /** function to fetch node buffer with other process; */
    int (*fetch_node_buffer)(struct dg_phys_info *phys, MPI_Request *send_requests, MPI_Request *recv_requests);
    /** function to fetch cell buffer with other process; */
    int (*fetch_cell_buffer)(struct dg_phys_info *phys, MPI_Request *send_requests, MPI_Request *recv_requests);

    void (*cell_mean)(struct dg_phys_info *phys);
}dg_phys_info;

/** slip wall and non-slop wall condition */
typedef int (*Wall_Condition)(dg_real nx, dg_real ny, dg_real *varM, dg_real *varP);
/** open boundary condition */
typedef int (*OBC_Fun)(dg_real nx, dg_real ny, dg_real *f_M, dg_real *f_ext, int obc_ind, dg_real *f_P);
/** two dimensional flux terms */
typedef int (*Nodal_Flux_Fun)(dg_real *var, dg_real *Eflux, dg_real *Gflux);
/** two dimensional numberical flux function */
typedef int (*Numerical_Flux)(dg_real nx, dg_real ny, dg_real *varM, dg_real *varP, dg_real *Fhs);

/** creator functions */
dg_phys_info* dg_phys_info_create(int Nfields, dg_area *area);
/** free function */
void dg_phys_info_free(dg_phys_info *phys_info);

#define dg_phys_info_area(info) info->area
#define dg_phys_info_cell(info) dg_area_cell(info->area)
#define dg_phys_info_grid(info) dg_area_grid(info->area)
#define dg_phys_info_region(info) dg_area_region(info->area)
#define dg_phys_info_mesh(info) dg_area_mesh(info->area)
#define dg_phys_info_edge(info) dg_area_edge(info->area)
#define dg_phys_info_procid(info) dg_area_procid(info->area)
#define dg_phys_info_nprocs(info) dg_area_nprocs(info->area)

#endif //DGOM_DG_PHYS_INFO_H