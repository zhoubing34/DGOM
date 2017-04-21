#ifndef DGOM_PHYSDOMAIN_H
#define DGOM_PHYSDOMAIN_H

#include "MultiArea/Edge/dg_edge.h"
#include "dg_phys_obc.h"
#include "dg_phys_info.h"
#include "dg_phys_LDG.h"
#include "PhysField/Limiter/dg_phys_limiter.h"


/**
 * @brief physical field structure.
 */
typedef struct dg_phys{
    dg_phys_info *info; ///< pointer to dg_phys_info structure;
    dg_phys_obc *obc; ///< pointer to dg_phys_obc structure;
    dg_phys_limiter *limiter; ///< pointer to dg_phys_limiter structure;
    dg_phys_LDG *ldg; ///< pointer to dg_phys_LDG structure;

    /** cell mean */
    void (*cell_mean)(struct dg_phys *phys);
    /** initialize from input file */
    void (*initialize_from_file)(struct dg_phys *phys, char *filename);
    /** function to fetch node buffer with other process */
    int (*fetch_node_buffer)(struct dg_phys *phys, MPI_Request *send_requests, MPI_Request *recv_requests);
    /** function to fetch cell buffer with other process */
    int (*fetch_cell_buffer)(struct dg_phys *phys, MPI_Request *send_requests, MPI_Request *recv_requests);
    /** add NetCDF format boundary condition file to physical field */
    void (*attach_obc_ncfile)(struct dg_phys *phys, char *filename);
    /** obtain the interpolated open boundary data from file */
    void (*update_obc_data)(struct dg_phys *phys, double elapseTime);

    void (*set_limiter)(struct dg_phys *phys, Limiter_Type type);
    void (*set_indicator)(struct dg_phys *phys, Indicator_Type type);
    /** use slope limiter to eliminate the oscillatory */
    void (*limit)(struct dg_phys *phys, double parameter);
} dg_phys;

dg_phys* dg_phys_create(int Nfields, dg_area *area);
void dg_phys_free(dg_phys *phys);


#define dg_phys_cell(phys) dg_phys_info_cell(phys->info)
#define dg_phys_grid(phys) dg_phys_info_grid(phys->info)
#define dg_phys_region(phys) dg_phys_info_region(phys->info)
#define dg_phys_mesh(phys) dg_phys_info_mesh(phys->info)
#define dg_phys_edge(phys) dg_phys_info_edge(phys->info)
#define dg_phys_nprocs(phys) dg_area_nprocs(dg_phys_info_area(phys->info))
#define dg_phys_procid(phys) dg_area_procid(dg_phys_info_area(phys->info))

#define dg_phys_Nfield(phys) (phys->info->Nfield)
#define dg_phys_f_Q(phys) (phys->info->f_Q)
#define dg_phys_f_recvQ(phys) (phys->info->f_recvQ)
#define dg_phys_c_Q(phys) phys->info->c_Q
#define dg_phys_c_recvQ(phys) phys->info->c_recvQ

#define dg_phys_f_rhsQ(phys) (phys->info->f_rhsQ)
#define dg_phys_f_resQ(phys) phys->info->f_resQ

#define dg_phys_f_extQ(phys) (phys->obc->f_extQ)

#endif //DGOM_PHYSDOMAIN_H
