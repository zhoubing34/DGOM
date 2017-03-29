//
// Created by li12242 on 17/3/29.
//
#include "swe_lib.h"

int swe_flux_term(dg_real *var, dg_real *Eflux, dg_real *Gflux){
    const dg_real h  = var[0];
    const dg_real qx = var[1];
    const dg_real qy = var[2];

    extern SWE_Solver solver;
    const dg_real hcrit = solver.hcrit;
    const dg_real gra = solver.gra;

    if(h>hcrit){
        Eflux[0] = qx; Eflux[1] = 0.5*gra*h*h + qx*qx/h; Eflux[2] = qx*qy/h;
        Gflux[0] = qy; Gflux[1] = qx*qy/h; Gflux[2] = qy*qy/h + 0.5*gra*h*h;
        Eflux[3] = 0;
        Gflux[3] = 0;
    }else{
        Eflux[0] = 0; Eflux[1] = 0; Eflux[2] = 0; Eflux[3] = 0;
        Gflux[0] = 0; Gflux[1] = 0; Gflux[2] = 0; Gflux[3] = 0;
    }
    return 0;
}

void swe_hll_cs(dg_real *varM, dg_real *varP, dg_real *sM, dg_real *sP){

    extern SWE_Solver solver;
    const dg_real hcrit = solver.hcrit;
    const dg_real gra = solver.gra;

    const dg_real hM = varM[0], hP = varP[0];
    const dg_real qnM = varM[1], qnP = varP[1];
    dg_real unM, unP;

    if( (hM>hcrit) & (hP>hcrit) ){ // both wet
        unM=qnM/hM;
        unP=qnP/hP;
        dg_real us = (dg_real)(0.5*(unM + unP)   + sqrt(gra*hM)   - sqrt(gra*hP));
        dg_real cs = (dg_real)(0.5*(sqrt(gra*hM) + sqrt(gra*hP) ) + 0.25*(unM - unP));

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

int swe_hll_flux(dg_real nx, dg_real ny, dg_real *varM, dg_real *varP, dg_real *Fhs){

    const dg_real hM  = varM[0], hP = varP[0];
    const dg_real qxM = varM[1], qxP = varP[1];
    const dg_real qyM = varM[2], qyP = varP[2];

    dg_real fn_M[3], fn_P[3];

    fn_M[0] = hM;
    fn_M[1] =  qxM*nx + qyM*ny; // local normal flux
    fn_M[2] = -qxM*ny + qyM*nx; // local tangential flux

    fn_P[0] = hP;
    fn_P[1] =  qxP*nx + qyP*ny; // next normal flux
    fn_P[2] = -qxP*ny + qyP*nx; // next tangential flux

    dg_real sM, sP;
    swe_hll_cs(fn_M, fn_P, &sM, &sP);

    dg_real eflux_M[3], eflux_P[3];
    dg_real gflux_M[3], gflux_P[3];

    swe_flux_term(fn_M, eflux_M, gflux_M);
    swe_flux_term(fn_P, eflux_P, gflux_P);

    dg_real Fhn=0, Fqxn=0, Fqyn=0;
    /* HLL function */
    if ( (sM>=0) & (sP>0) ){
        Fhn = eflux_M[0]; Fqxn = eflux_M[1]; Fqyn = eflux_M[2];
    }else if((sM<0) & (sP>0)){
        Fhn  = (sP*eflux_M[0] - sM*eflux_P[0]  + sM*sP*(hP  - hM ))/(sP - sM);
        Fqxn = (sP*eflux_M[1] - sM*eflux_P[1] + sM*sP*(fn_P[1] - fn_M[1]))/(sP - sM);
        Fqyn = (sP*eflux_M[2] - sM*eflux_P[2] + sM*sP*(fn_P[2] - fn_M[2]))/(sP - sM);
    }else if( (sM<0)&(sP<=0) ){
        Fhn = eflux_P[0]; Fqxn = eflux_P[1]; Fqyn = eflux_P[2];
    }else if( fabs(sM-0)<EPS & fabs(sP-0)<EPS ){
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
    Fhs[3] = 0;
    return 0;
}