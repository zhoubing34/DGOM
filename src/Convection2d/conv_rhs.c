#include <MultiRegions/Mesh/dg_mesh.h>
#include "conv_driver2d.h"
#include "PhysField/pf_fetchBuffer.h"
#include "PhysField/pf_strong_volume_flux2d.h"
#include "PhysField/pf_strong_surface_flux2d.h"
#include "PhysField/pf_strong_viscosity_LDG_flux2d.h"

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


void conv_rhs(physField *phys, dg_real frka, dg_real frkb, dg_real fdt){

    dg_mesh *mesh = phys->mesh;
    const int K = mesh->grid->K;
    const int nprocs = mesh->nprocs;
    const int Np = phys->cell->Np;
    const int Nfield = phys->Nfield;

    extern conv_solver2d solver;

    /* mpi request buffer */
    MPI_Request *mpi_send_requests = (MPI_Request *) calloc((size_t) nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_recv_requests = (MPI_Request *) calloc((size_t) nprocs, sizeof(MPI_Request));

    int Nmess=0;
    /* fetch nodal value through all procss */
    pf_fetchNodeBuffer2d(phys, mpi_send_requests, mpi_recv_requests, &Nmess);

    /* volume integral */
    pf_strong_volume_flux2d(phys, conv_fluxTerm);

    /* waite to recv */
    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_recv_requests, instatus);

    /* surface integral */
    pf_strong_surface_flux2d(phys, NULL, NULL, conv_fluxTerm, conv_upWindFlux);

    /* waite for finishing send buffer */
    MPI_Waitall(Nmess, mpi_send_requests, instatus);

    /* viscosity flux */
    if(solver.viscosity > DIFF_THRESHOLD) {
        dg_real c11 = (dg_real) solver.LDG_parameter[0];
        dg_real c12 = (dg_real) solver.LDG_parameter[1];
        dg_real c22 = (dg_real) solver.LDG_parameter[2];
        pf_strong_viscosity_LDG_flux2d(phys, NULL, NULL, c11, c12, c22);
    }

    dg_real *f_resQ = phys->f_resQ;
    dg_real *f_Q = phys->f_Q;
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