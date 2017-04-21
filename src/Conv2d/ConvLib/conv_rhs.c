
#include "conv_lib2d.h"

static int conv_flux(dg_real *var, dg_real *Eflux, dg_real *Gflux){
    const dg_real c = var[0];
    const dg_real u = var[1];
    const dg_real v = var[2];

    Eflux[0] = u*c; Eflux[1] = 0; Eflux[2] = 0;
    Gflux[0] = v*c; Gflux[1] = 0; Gflux[2] = 0;

    return 0;
}

static int conv_upwind_flux(dg_real nx, dg_real ny, dg_real *f_M, dg_real *f_P, dg_real *Fhs){
    const dg_real cM = f_M[0], cP = f_P[0];
    const dg_real uM = f_M[1], uP = f_P[1];
    const dg_real vM = f_M[2], vP = f_P[2];

    if( (uM*nx+vM*ny) > 0.0 ){
        Fhs[0] = nx*cM*uM + ny*cM*vM;
    }else{
        Fhs[0] = nx*cP*uP + ny*cP*vP;
    }

    Fhs[1] = 0; Fhs[2] = 0;
    return 0;
}

static int conv_obc_func(dg_real nx, dg_real ny, dg_real *f_M, dg_real *f_ext, int obc_ind, dg_real *f_P){
    const dg_real cM = f_M[0];
    const dg_real uM = f_M[1];
    const dg_real vM = f_M[2];

    if( (uM*nx+vM*ny) > 0.0 ){
        f_P[0] = cM;
    }else{
        f_P[0] = 0 - cM;
    }

    f_P[1] = 0; f_P[2] = 0;
    return 0;
}

static int vis_func(dg_real *f_Q, dg_real *vis){
    vis[0] = f_Q[0];
    vis[1] = 0;
    vis[2] = 0;
    return 0;
}


void conv_rhs(dg_phys *phys, dg_real frka, dg_real frkb, dg_real fdt){
    dg_grid *grid = dg_phys_grid(phys);
    const int K = dg_grid_K(grid);
    const int nprocs = dg_phys_nprocs(phys);
    const int Np = dg_cell_Np(dg_phys_cell(phys));
    const int Nfield = dg_phys_Nfield(phys);

    dg_real *f_Q = dg_phys_f_Q(phys);

    /* mpi request buffer */
    MPI_Request *mpi_send_requests = (MPI_Request *) calloc((size_t) nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_recv_requests = (MPI_Request *) calloc((size_t) nprocs, sizeof(MPI_Request));

    /* fetch nodal value through all procss */
    int Nmess = phys->fetch_node_buffer(phys, mpi_send_requests, mpi_recv_requests);

    /* volume vol_integral */
    dg_phys_strong_vol_opt2d(phys, conv_flux);

    /* waite to recv */
    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_recv_requests, instatus);

    /* surface vol_integral */
    dg_phys_strong_surf_opt2d(phys, NULL, NULL, conv_obc_func, conv_flux, conv_upwind_flux);

    /* waite for finishing send buffer */
    MPI_Waitall(Nmess, mpi_send_requests, instatus);

    /* viscosity flux */
    extern Conv_Solver solver;
    if(solver.vis_flag){
        dg_phys_LDG_solve_vis_opt2d(phys, vis_func, NULL, NULL, conv_obc_func);
    }

    dg_real *f_resQ = dg_phys_f_resQ(phys);
    const dg_real *f_rhsQ = dg_phys_f_rhsQ(phys);
    int t;
    for(t=0;t<K*Np*Nfield;++t){
        f_resQ[t] = frka*f_resQ[t]+fdt*f_rhsQ[t];   // calculate the resdiual of the equation
        f_Q[t]   += frkb*f_resQ[t];                 // evaluate scalar at next internal time step
    }

    /* make sure all messages went out */
    MPI_Status outstatus[nprocs];
    MPI_Waitall(Nmess, mpi_send_requests, outstatus);

    free(mpi_recv_requests);
    free(mpi_send_requests);
}