//
// Created by li12242 on 17/1/14.
//

#include <PhysField/pf_phys.h>
#include "swe_rhs.h"
#include "PhysField/pf_fetchBuffer.h"
#include "PhysField/pf_strong_volume_flux2d.h"
#include "PhysField/pf_strong_surface_flux2d.h"
#include "swe_driver2d.h"

#define TOL 1e-8
static real gra;
static real hcrit;

static void swe_sour(swe_solver *solver, physField *phys);

int swe_flux_term(real *var, real *Eflux, real *Gflux){
    const real h  = var[0];
    const real qx = var[1];
    const real qy = var[2];

    if(h>hcrit){
        Eflux[0] = qx; Eflux[1] = 0.5*gra*h*h + qx*qx/h; Eflux[2] = qx*qy/h;
        Gflux[0] = qy; Gflux[1] = qx*qy/h; Gflux[2] = qy*qy/h + 0.5*gra*h*h;
    }else{
        Eflux[0] = 0; Eflux[1] = 0; Eflux[2] = 0;
        Gflux[0] = 0; Gflux[1] = 0; Gflux[2] = 0;
    }
    return 0;
}


void swe_hll_cs(real *varM, real *varP, real *sM, real *sP){

    const real hM = varM[0], hP = varP[0];
    const real qnM = varM[1], qnP = varP[1];
    real unM, unP;

    if( (hM>hcrit) & (hP>hcrit) ){ // both wet
        unM=qnM/hM;
        unP=qnP/hP;
        real us = (real)(0.5*(unM + unP)   + sqrt(gra*hM)   - sqrt(gra*hP));
        real cs = (real)(0.5*(sqrt(gra*hM) + sqrt(gra*hP) ) + 0.25*(unM - unP));

        *sM = min(unM-sqrt(gra*hM), us-cs);
        *sP = max(unP+sqrt(gra*hP), us+cs);
    }else if ( (hM>hcrit) & (hP<=hcrit) ){
        unM=qnM/hM;
        *sM = (unM -  sqrt(gra*hM) );
        *sP = (unM +2*sqrt(gra*hM) );
    }else if ( (hM<=hcrit) & (hP>hcrit) ){
        unP=qnP/hP;
        *sM = (unP -2*sqrt(gra*hP) );
        *sP = (unP +  sqrt(gra*hP) );
    }else{ /* both dry element */
        *sM = 0; *sP = 0;
    }
    return;
}

int swe_hll_flux(real nx, real ny, real *varM, real *varP, real *Fhs){

    const real hM  = varM[0], hP = varP[0];
    const real qxM = varM[1], qxP = varP[1];
    const real qyM = varM[2], qyP = varP[2];

    real fn_M[3], fn_P[3];

    fn_M[0] = hM;
    fn_M[1] =  qxM*nx + qyM*ny; // local normal flux
    fn_M[2] = -qxM*ny + qyM*nx; // local tangential flux

    fn_P[0] = hP;
    fn_P[1] =  qxP*nx + qyP*ny; // next normal flux
    fn_P[2] = -qxP*ny + qyP*nx; // next tangential flux

    real sM, sP;
    swe_hll_cs(fn_M, fn_P, &sM, &sP);

    real eflux_M[3], eflux_P[3];
    real gflux_M[3], gflux_P[3];

    swe_flux_term(fn_M, eflux_M, gflux_M);
    swe_flux_term(fn_P, eflux_P, gflux_P);

    real Fhn=0, Fqxn=0, Fqyn=0;
    /* HLL function */
    if ( (sM>=0) & (sP>0) ){
        Fhn = eflux_M[0]; Fqxn = eflux_M[1]; Fqyn = eflux_M[2];
    }else if((sM<0) & (sP>0)){
        Fhn  = (sP*eflux_M[0] - sM*eflux_P[0]  + sM*sP*(hP  - hM ))/(sP - sM);
        Fqxn = (sP*eflux_M[1] - sM*eflux_P[1] + sM*sP*(fn_P[1] - fn_M[1]))/(sP - sM);
        Fqyn = (sP*eflux_M[2] - sM*eflux_P[2] + sM*sP*(fn_P[2] - fn_M[2]))/(sP - sM);
    }else if( (sM<0)&(sP<=0) ){
        Fhn = eflux_P[0]; Fqxn = eflux_P[1]; Fqyn = eflux_P[2];
    }else if( fabs(sM-0)<TOL & fabs(sP-0)<TOL ){
        Fhn = 0; Fqxn = 0; Fqyn = 0;
    }else{
        fprintf(stderr, "swe_hll_flux (%s): %d\n"
                        "flux speed error sM=%e, sP=%e\n"
                        "varM=[%e, %e, %e]\n"
                        "varP=[%e, %e, %e]\n",
                __FILE__, __LINE__, sM, sP,
                varM[0],varM[1],varM[2],
                varP[0],varP[1],varP[2]);
//        MPI_Abort(MPI_COMM_WORLD, -1);
        return 1;
    }
    Fhs[0] = Fhn;
    Fhs[1] = Fqxn*nx - Fqyn*ny;
    Fhs[2] = Fqxn*ny + Fqyn*nx;
    return 0;
}

int swe_slip_wall(real nx, real ny, real *varM, real *varP){
    const real hM  = varM[0];
    const real qxM = varM[1];
    const real qyM = varM[2];
    real qnM =  qxM*nx + qyM*ny; // outward normal flux
    real qvM = -qxM*ny + qyM*nx; // outward tangential flux
    // adjacent value
    varP[0] = hM;
    varP[1] = (-qnM)*nx - qvM*ny;
    varP[2] = (-qnM)*ny + qvM*nx;
    return 0;
}

int swe_nonslip_wall(real nx, real ny, real *varM, real *varP){
    const real hM  = varM[0];
    const real qxM = varM[1];
    const real qyM = varM[2];
    // adjacent value
    varP[0] = hM;
    varP[1] = -qxM;
    varP[2] = -qyM;
    return 0;
}

void swe_rhs(swe_solver *solver, real frka, real frkb, real fdt){

    physField *phys = solver->phys;
    const int K = phys->grid->K;
    const int Np = phys->cell->Np;
    const int Nfield = phys->Nfield;
    const int nprocs = phys->mesh->nprocs;

    extern real gra, hcrit;
    gra = solver->gra;
    hcrit = solver->hcrit;

    /* mpi request buffer */
    MPI_Request *mpi_send_requests = (MPI_Request *) calloc((size_t) nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_recv_requests = (MPI_Request *) calloc((size_t) nprocs, sizeof(MPI_Request));

    int Nmess=0;
    /* fetch nodal value through all procss */
    pf_fetchNodeBuffer2d(phys, mpi_send_requests, mpi_recv_requests, &Nmess);

    /* volume integral */
    pf_strong_volume_flux2d(phys, swe_flux_term);

    /* waite to recv */
    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_recv_requests, instatus);

    /* surface integral */
    pf_strong_surface_flux2d(phys, swe_slip_wall, swe_nonslip_wall, swe_flux_term, swe_hll_flux);

    /* waite for finishing send buffer */
    MPI_Waitall(Nmess, mpi_send_requests, instatus);

    /* viscosity flux */
    /* source term */
    swe_sour(solver, phys);

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

    free(mpi_recv_requests);
    free(mpi_send_requests);

    return;
}

/**
 * @brief add source term to the rhs field
 */
static void swe_sour(swe_solver *solver, physField *phys){

    const real gra     = (real)solver->gra;
    const real *f_Dr   = phys->cell->f_Dr;
    const real *f_Ds   = phys->cell->f_Ds;
    const real *vgeo = phys->vgeo;
    const real hcrit = solver->hcrit;
    const int K        = phys->grid->K;
    const int Np       = phys->cell->Np;
    const int Nfields  = phys->Nfield;

    const real M2 = (real) solver->roughness*solver->roughness;

    real *f_Q = phys->f_Q;
    real *f_rhsQ = phys->f_rhsQ;

    const double p = 10/3;

    register int k,n,m,geoid=0;

    for(k=0;k<K;k++){
        real *bot = solver->bot + k*Np;
        for (n=0;n<Np;++n) {

            const real *ptDr = f_Dr + n*Np;
            const real *ptDs = f_Ds + n*Np;
            const real drdx  = vgeo[geoid++], drdy = vgeo[geoid++];
            const real dsdx  = vgeo[geoid++], dsdy = vgeo[geoid++];

            const int ind = (k*Np+n)*Nfields;
            const real h = f_Q[ind];
            const real qx = f_Q[ind];
            const real qy = f_Q[ind];
            const real qn = sqrt(qx*qx+qy*qy);

            real rhs_qx = 0, rhs_qy = 0;
            if(h>hcrit){ // for wet nodes
                for (m=0;m<Np;++m) {
                    const real dr = ptDr[m];
                    const real ds = ptDs[m];
                    const real dx = drdx*dr + dsdx*ds;
                    const real dy = drdy*dr + dsdy*ds;

                    const real b = (real) bot[m];

                    rhs_qx += dx * b;
                    rhs_qy += dy * b;
                }

                rhs_qx += qx*qn*M2/pow(h, p);
                rhs_qy += qy*qn*M2/pow(h, p);
            }

            f_rhsQ[ind+1] += -gra*h*rhs_qx;
            f_rhsQ[ind+2] += -gra*h*rhs_qy;
        }
    }
    return;
}