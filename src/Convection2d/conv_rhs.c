
#include "conv_driver.h"

int conv_fluxTerm(dg_real *var, dg_real *Eflux, dg_real *Gflux){
    const dg_real c = var[0];
    const dg_real u = var[1];
    const dg_real v = var[2];

    Eflux[0] = u*c; Eflux[1] = 0; Eflux[2] = 0;
    Gflux[0] = v*c; Gflux[1] = 0; Gflux[2] = 0;

    return 0;
}

int conv_upWindFlux(dg_real nx, dg_real ny, dg_real *varM, dg_real *varP, dg_real *Fhs){
    const dg_real cM = varM[0], cP = varP[0];
    const dg_real uM = varM[1], uP = varP[1];
    const dg_real vM = varM[2], vP = varP[2];

    if( (uM*nx+vM*ny) > 0.0 ){
        Fhs[0] = nx*cM*uM + ny*cM*vM;
    }else{
        Fhs[0] = nx*cP*uP + ny*cP*vP;
    }

    Fhs[1] = 0; Fhs[2] = 0;
    return 0;
}


void conv_rhs(dg_phys *phys, dg_real frka, dg_real frkb, dg_real fdt){
    const int K = dg_grid_K(phys->grid);
    const int nprocs = dg_grid_nprocs(phys->grid);
    const int Np = dg_cell_Np(phys->cell);
    const int Nfield = dg_phys_Nfield(phys);

    dg_real *f_Q = phys->f_Q;

    /* mpi request buffer */
    MPI_Request *mpi_send_requests = (MPI_Request *) calloc((size_t) nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_recv_requests = (MPI_Request *) calloc((size_t) nprocs, sizeof(MPI_Request));

    /* fetch nodal value through all procss */
    int Nmess = phys->fetch_node_buffer(phys, mpi_send_requests, mpi_recv_requests);

    /* volume integral */
    dg_phys_strong_vol_opt2d(phys, conv_fluxTerm);

    /* waite to recv */
    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_recv_requests, instatus);

    /* surface integral */
    dg_phys_strong_surf_opt2d(phys, NULL, NULL, conv_fluxTerm, conv_upWindFlux);

    /* waite for finishing send buffer */
    MPI_Waitall(Nmess, mpi_send_requests, instatus);

    /* viscosity flux */

    dg_real *f_resQ = phys->f_resQ;
    const dg_real *f_rhsQ = phys->f_rhsQ;
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