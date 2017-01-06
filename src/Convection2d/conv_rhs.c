#include <MultiRegions/mr_mesh.h>
#include "conv_driver2d.h"
#include "PhysField/phys_fetchBuffer.h"
#include "PhysField/phys_strong_volume_flux2d.h"
#include "PhysField/phys_strong_surface_flux2d.h"
#include "PhysField/phys_strong_viscosity_LDG_flux2d.h"

int conv_fluxTerm(real *var, real *Eflux, real *Gflux){
    const real c = var[0];
    const real u = var[1];
    const real v = var[2];

    Eflux[0] = u*c; Eflux[1] = 0; Eflux[2] = 0;
    Gflux[0] = v*c; Gflux[1] = 0; Gflux[2] = 0;

    return 0;
}

int conv_upWindFlux(real nx, real ny, real *varM, real *varP, real *Fhs){
    const real cM = varM[0], cP = varP[0];
    const real uM = varM[1], uP = varP[1];
    const real vM = varM[2], vP = varP[2];

    if( (uM*nx+vM*ny) > 0.0 ){
        Fhs[0] = nx*cM*uM + ny*cM*vM;
    }else{
        Fhs[0] = nx*cP*uP + ny*cP*vP;
    }

    Fhs[1] = 0; Fhs[2] = 0;
    return 0;
}


void conv_rhs(physField *phys, float frka, float frkb, float fdt){

    parallMesh *mesh = phys->mesh;
    const int K = mesh->grid->K;
    const int nprocs = mesh->nprocs;
    const int Np = phys->cell->Np;
    const int Nfield = phys->Nfield;

    extern conv_solver2d solver;

    /* mpi request buffer */
    MPI_Request mpi_send_requests[nprocs];
    MPI_Request mpi_recv_requests[nprocs];

    int Nmess=0;
    /* fetch nodal value through all procss */
    phys_fetchNodeBuffer2d(phys, mpi_send_requests, mpi_recv_requests, &Nmess);

    /* volume integral */
    phys_strong_volume_flux2d(phys, conv_fluxTerm);

    /* waite to recv */
    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_recv_requests, instatus);

    /* surface integral */
    phys_strong_surface_integral2d(phys, NULL, NULL, conv_fluxTerm, conv_upWindFlux);

    /* waite for finishing send buffer */
    MPI_Waitall(Nmess, mpi_send_requests, instatus);

    /* viscosity flux */
    if(solver.caseid == conv_advection_diffusion)
        phys_strong_viscosity_LDG_flux2d(phys, NULL, NULL, 0, 0, 0);


    real *f_resQ = phys->f_resQ;
    real *f_Q = phys->f_Q;
    const real *f_rhsQ = phys->f_rhsQ;
    int t;
    for(t=0;t<K*Np*Nfield;++t){
        f_resQ[t] = frka*f_resQ[t]+fdt*f_rhsQ[t];   // calculate the resdiual of the equation
        f_Q[t]   += frkb*f_resQ[t];                 // evaluate scalar at next internal time step
    }

    /* make sure all messages went out */
    MPI_Status outstatus[nprocs];
    MPI_Waitall(Nmess, mpi_send_requests, outstatus);
}