//
// Created by li12242 on 17/1/1.
//

#include <StandCell/dg_cell.h>
#include <MultiRegions/Grid/dg_grid.h>
#include "pf_strong_viscosity_LDG_flux2d.h"
#include "pf_fetchBuffer.h"
#include "pf_phys.h"

#define DEBUG 0

#if DEBUG
#include "LibUtilities/UTest.h"
#endif

static void phys_viscosityflux2d(physField *phys,
                                 viscosity_wall_condition_func slipwall_condition,
                                 viscosity_wall_condition_func non_slipwall_condition,
                                 dg_real c11, dg_real c12, dg_real c22);

static void phys_auxiliaryflux2d(physField *phys,
                                 viscosity_wall_condition_func slipwall_condition,
                                 viscosity_wall_condition_func non_slipwall_condition,
                                 dg_real c11, dg_real c12, dg_real c22);

void pf_strong_viscosity_LDG_flux2d(physField *phys,
                                    viscosity_wall_condition_func slipwall_condition,
                                    viscosity_wall_condition_func non_slipwall_condition,
                                    dg_real c11, dg_real c12, dg_real c22){

    // check LDG solver is initialized
    if(phys->viscosity == NULL){
        fprintf(stderr, "PhysField (%s): line %d\n"
                "The LDG solver for viscosity is not allocated\n", __FILE__, __LINE__);
    }

    const int procid = phys->grid->procid;
    const int nprocs = phys->grid->nprocs;
    const int Nfp = phys->cell->Nfp;
    const int Nfield = phys->Nfield;

    dg_real *p_Q = phys->viscosity->px_Q;
    dg_real *q_Q = phys->viscosity->py_Q;
    dg_real *p_inQ = phys->viscosity->px_inQ;
    dg_real *q_inQ = phys->viscosity->py_inQ;
    dg_real *p_outQ = phys->viscosity->px_outQ;
    dg_real *q_outQ = phys->viscosity->py_outQ;

    register int n;

    // calculate the auxiliary variable
    phys_auxiliaryflux2d(phys, slipwall_condition, non_slipwall_condition, c11, c12, c22);



    // fetch p_Q and q_Q buffer
    for(n=0;n<phys->parallNodeNum;++n){
        p_outQ[n] = p_Q[phys->nodeIndexOut[n]];
        q_outQ[n] = q_Q[phys->nodeIndexOut[n]];
    }

    int Nout[nprocs];
    for(n=0;n<nprocs;n++){
        Nout[n] = phys->mesh->Parf[n]*Nfield*Nfp;
    }

    MPI_Request mpi_send_requests[nprocs], mpi_recv_requests[nprocs];
    int Nmess=0;

    /* do sends and recv for p_Q */
    pf_fetchBuffer(procid, nprocs, Nout, p_outQ, p_inQ,
                   mpi_send_requests, mpi_recv_requests, &Nmess);

    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_recv_requests, instatus);
    MPI_Waitall(Nmess, mpi_send_requests, instatus);

    /* do sends and recv for q_Q */
    pf_fetchBuffer(procid, nprocs, Nout, q_outQ, q_inQ,
                   mpi_send_requests, mpi_recv_requests, &Nmess);
    MPI_Waitall(Nmess, mpi_recv_requests, instatus);
    MPI_Waitall(Nmess, mpi_send_requests, instatus);

    // calculate the RHS of physical field
    phys_viscosityflux2d(phys, slipwall_condition, non_slipwall_condition, c11, c12, c22);

}


static void phys_auxiliaryflux2d(physField *phys,
                                 viscosity_wall_condition_func slipwall_condition,
                                 viscosity_wall_condition_func non_slipwall_condition,
                                 dg_real c11, dg_real c12, dg_real c22){

    const int K = phys->grid->K;
    const int Np = phys->cell->Np;
    const int Nfp = phys->cell->Nfp;
    const int Nfaces = phys->cell->Nfaces;
    const int Nfield = phys->Nfield;

    dg_real *p_Q = phys->viscosity->px_Q;
    dg_real *q_Q = phys->viscosity->py_Q;
    dg_real *p_inQ = phys->viscosity->px_inQ;
    dg_real *q_inQ = phys->viscosity->py_inQ;

    dg_real *f_inQ = phys->f_inQ;
    dg_real *f_ext = phys->f_ext;

    dg_real *f_Q = phys->f_Q;
    dg_real *f_Dr = phys->cell->f_Dr;
    dg_real *f_Ds = phys->cell->f_Ds;
    dg_real *f_LIFT = phys->cell->f_LIFT;
    dg_real *vgeo = phys->vgeo;
    dg_real *surfinfo = phys->surfinfo;

    dg_real *vis_Q = phys->viscosity->vis_Q;

    register int m,n,k,fld,geoid=0,surfid=0;

    dg_real f_varM[Nfield], f_varP[Nfield];
    dg_real p_varM[Nfield], p_varP[Nfield], p_dflux[Nfp*Nfaces*Nfield];
    dg_real q_varM[Nfield], q_varP[Nfield], q_dflux[Nfp*Nfaces*Nfield];
    dg_real vis_varM[Nfield];

#if DEBUG
    FILE *fp = CreateLog("phys_auxiliaryflux2d_log",
                         phys->mesh->procid, phys->mesh->nprocs);
#endif

    for(k=0;k<K;k++){
        // volume integral for p and q
        for(n=0;n<Np;n++){
            const dg_real *ptDr = f_Dr+n*Np; // n-th row of Dr
            const dg_real *ptDs = f_Ds+n*Np; // n-th row of Ds

            const dg_real drdx = vgeo[geoid++]; // volume geometry for n-th point
            const dg_real drdy = vgeo[geoid++]; // volume geometry for n-th point
            const dg_real dsdx = vgeo[geoid++]; // volume geometry for n-th point
            const dg_real dsdy = vgeo[geoid++]; // volume geometry for n-th point
            const dg_real *visc = phys->viscosity->vis_Q+(k*Np+n)*Nfield;

            dg_real *p = p_Q + (n+k*Np)*Nfield;
            dg_real *q = q_Q + (n+k*Np)*Nfield;

            // initialization: set px and py to zero
            for(fld=0;fld<Nfield;fld++){
                p[fld] = 0;
                q[fld] = 0;
            }

            for(m=0;m<Np;++m){
                const dg_real dr = ptDr[m]; // m-th column for Dr
                const dg_real ds = ptDs[m]; // m-th column for Ds
                const dg_real dx = drdx*dr+dsdx*ds;
                const dg_real dy = drdy*dr+dsdy*ds;

                const dg_real *c = f_Q + (m+k*Np)*Nfield;

                for(fld=0;fld<Nfield;fld++){
                    p[fld] += visc[fld]*dx*c[fld];
                    q[fld] += visc[fld]*dy*c[fld];
                }
            }
        }

        // surface integral
        for(n=0;n<Nfp*Nfaces;n++){
            int  idM = (int)surfinfo[surfid++];
            int  idP = (int)surfinfo[surfid++];
            const dg_real fsc = surfinfo[surfid++];
            const int  bstype = (int)surfinfo[surfid++];
            const dg_real nx = surfinfo[surfid++];
            const dg_real ny = surfinfo[surfid++];

            // local face2d values
            for(fld=0;fld<Nfield;fld++){
                f_varM[fld] = f_Q[idM];
                p_varM[fld] = p_Q[idM];
                q_varM[fld] = q_Q[idM];
                vis_varM[fld] = vis_Q[idM++];
            }

#if DEBUG
            fprintf(fp, "\nNfp=%d, bctype=%d, ", n, bstype);
#endif

            // adjacent nodal values
            switch (bstype){
                case INNERLOC:
                    for(fld=0;fld<Nfield;fld++){
                        f_varP[fld] = f_Q[idP];
                        p_varP[fld] = p_Q[idP];
                        q_varP[fld] = q_Q[idP++];
                    }
                    break;
                case INNERBS:
                    for(fld=0;fld<Nfield;fld++){
                        f_varP[fld] = f_inQ[idP];
                        p_varP[fld] = p_inQ[idP];
                        q_varP[fld] = q_inQ[idP++];
                    }
                    break;
                case SLIPWALL:
                    slipwall_condition(nx, ny, f_varM, f_varP,
                                       p_varM, p_varP, q_varM, q_varP);
                    break;
                case NSLIPWALL:
                    non_slipwall_condition(nx, ny, f_varM, f_varP,
                                           p_varM, p_varP, q_varM, q_varP);
                    break;
                default: // open boundary condition
                    for(fld=0;fld<Nfield;fld++){
                        f_varP[fld] = f_ext[idP++];
                        p_varP[fld] = p_varM[fld]; // zero gradient
                        q_varP[fld] = q_varM[fld];
                    }
                    break;
            }

            dg_real *flux_P = p_dflux + n*Nfield;
            dg_real *flux_Q = q_dflux + n*Nfield;

            for(fld=0;fld<Nfield;fld++){
                const dg_real df = (f_varM[fld] - f_varP[fld]);
                const dg_real dp = nx*(p_varM[fld] - p_varP[fld]) + ny*(q_varM[fld] - q_varP[fld]);
#if DEBUG
                fprintf(fp, "%f, ", df);
#endif
                flux_P[fld] = vis_varM[fld]*fsc*( nx*( -df*0.5 + c12*nx*df - c22*dp ) );
                flux_Q[fld] = vis_varM[fld]*fsc*( ny*( -df*0.5 + c12*ny*df - c22*dp ) );
            }
        }

#if DEBUG
        PrintVector2File(fp, "elemental flux_P", p_dflux, Nfield*Nfp*Nfaces);
        PrintVector2File(fp, "elemental flux_Q", q_dflux, Nfield*Nfp*Nfaces);
#endif

        for(n=0;n<Np;n++){
            const dg_real *ptLIFT = f_LIFT + n*Nfp*Nfaces;

            dg_real *p_rhsQ = p_Q + Nfield*(n+k*Np);
            dg_real *q_rhsQ = q_Q + Nfield*(n+k*Np);
            for(m=0;m<Nfp*Nfaces;m++){
                const dg_real L = ptLIFT[m];
                dg_real *flux_P = p_dflux+m*Nfield;
                dg_real *flux_Q = q_dflux+m*Nfield;

                for(fld=0;fld<Nfield;fld++){
                    p_rhsQ[fld] += L*flux_P[fld];
                    q_rhsQ[fld] += L*flux_Q[fld];
                }

            }
        }
    }
#if DEBUG
    fclose(fp);
#endif

}

static void phys_viscosityflux2d(physField *phys,
                                 viscosity_wall_condition_func slipwall_condition,
                                 viscosity_wall_condition_func non_slipwall_condition,
                                 dg_real c11, dg_real c12, dg_real c22){
    const int K = phys->grid->K;
    const int Np = phys->cell->Np;
    const int Nfp = phys->cell->Nfp;
    const int Nfaces = phys->cell->Nfaces;
    const int Nfield = phys->Nfield;

    dg_real *f_Q = phys->f_Q;
    dg_real *f_Dr = phys->cell->f_Dr;
    dg_real *f_Ds = phys->cell->f_Ds;
    dg_real *f_LIFT = phys->cell->f_LIFT;
    dg_real *f_rhs = phys->f_rhsQ;
    dg_real *vgeo = phys->vgeo;
    dg_real *surfinfo = phys->surfinfo;

    dg_real *p_Q = phys->viscosity->px_Q;
    dg_real *q_Q = phys->viscosity->py_Q;
    dg_real *p_inQ = phys->viscosity->px_inQ;
    dg_real *q_inQ = phys->viscosity->py_inQ;

    dg_real *f_inQ = phys->f_inQ;
    dg_real *f_ext = phys->f_ext;

    register int k,n,m,fld,geoid=0,surfid=0;

    dg_real f_varM[Nfield], f_varP[Nfield], f_dflux[Nfp*Nfaces*Nfield];
    dg_real p_varM[Nfield], p_varP[Nfield];
    dg_real q_varM[Nfield], q_varP[Nfield];

#if DEBUG
    FILE *fp = CreateLog("phys_viscosityflux2d_log",
                         phys->mesh->procid, phys->mesh->nprocs);
#endif

    for(k=0;k<K;k++){
        // volume integral for p and q
        for(n=0;n<Np;n++){
            const dg_real *ptDr = f_Dr+n*Np; // n-th row of Dr
            const dg_real *ptDs = f_Ds+n*Np; // n-th row of Ds

            const dg_real drdx = vgeo[geoid++]; // volume geometry for n-th point
            const dg_real drdy = vgeo[geoid++]; // volume geometry for n-th point
            const dg_real dsdx = vgeo[geoid++]; // volume geometry for n-th point
            const dg_real dsdy = vgeo[geoid++]; // volume geometry for n-th point

            dg_real *rhs = f_rhs + (k*Np+n)*Nfield;

            for(m=0;m<Np;++m){
                const dg_real dr = ptDr[m]; // m-th column for Dr
                const dg_real ds = ptDs[m]; // m-th column for Ds
                const dg_real dx = drdx*dr+dsdx*ds;
                const dg_real dy = drdy*dr+dsdy*ds;

                const dg_real *p = p_Q + (m+k*Np)*Nfield;
                const dg_real *q = q_Q + (m+k*Np)*Nfield;

                for(fld=0;fld<Nfield;fld++){
                    rhs[fld] += dx*p[fld] + dy*q[fld];
                }
            }
        }

        for(n=0;n<Nfaces*Nfp;n++){
            int  idM = (int)surfinfo[surfid++];
            int  idP = (int)surfinfo[surfid++];
            const dg_real fsc = surfinfo[surfid++];
            const int  bstype = (int)surfinfo[surfid++];
            const dg_real nx = surfinfo[surfid++];
            const dg_real ny = surfinfo[surfid++];

            // local face2d values
            for(fld=0;fld<Nfield;fld++){
                f_varM[fld] = f_Q[idM];
                p_varM[fld] = p_Q[idM];
                q_varM[fld] = q_Q[idM++];
            }

#if DEBUG
            fprintf(fp, "\nNfp=%d, bctype=%d, ", n, bstype);
#endif

            // adjacent nodal values
            switch (bstype){
                case INNERLOC:
                    for(fld=0;fld<Nfield;fld++){
                        f_varP[fld] = f_Q[idP];
                        p_varP[fld] = p_Q[idP];
                        q_varP[fld] = q_Q[idP++];
                    }
                    break;
                case INNERBS:
                    for(fld=0;fld<Nfield;fld++){
                        f_varP[fld] = f_inQ[idP];
                        p_varP[fld] = p_inQ[idP];
                        q_varP[fld] = q_inQ[idP++];
                    }
                    break;
                case SLIPWALL:
                    slipwall_condition(nx, ny, f_varM, f_varP,
                                       p_varM, p_varP, q_varM, q_varP);
                    break;
                case NSLIPWALL:
                    non_slipwall_condition(nx, ny, f_varM, f_varP,
                                           p_varM, p_varP, q_varM, q_varP);
                    break;
                default: // open boundary condition
                    for(fld=0;fld<Nfield;fld++){
                        f_varP[fld] = f_ext[idP++];
                        p_varP[fld] = p_varM[fld]; // zero gradient
                        q_varP[fld] = q_varM[fld];
                    }
                    break;
            }

            dg_real *flux = f_dflux + n*Nfield;

            for(fld=0;fld<Nfield;fld++){
                const dg_real df = (f_varM[fld] - f_varP[fld]);
                const dg_real dp = (p_varM[fld] - p_varP[fld]);
                const dg_real dq = (q_varM[fld] - q_varP[fld]);
                const dg_real dpn = nx*dp + ny*dq;

#if DEBUG
                fprintf(fp, "%f, %f, %f, ", df, dp, dq);
#endif
                flux[fld] = fsc*( -nx*dp*0.5 - ny*dq*0.5 - c11*df - c12*(nx+ny)*dpn);
            }
        }

#if DEBUG
        PrintVector2File(fp, "elemental flux", f_dflux, Nfield*Nfp*Nfaces);
#endif

        for(n=0;n<Np;n++){
            const dg_real *ptLIFT = f_LIFT + n*Nfp*Nfaces;
            for(m=0;m<Nfp*Nfaces;m++){
                const dg_real L = ptLIFT[m];
                dg_real *flux_Q = f_dflux+m*Nfield;
                dg_real *f_rhsQ = phys->f_rhsQ + Nfield*(n+k*Np);

                for(fld=0;fld<Nfield;fld++)
                    f_rhsQ[fld] += L*flux_Q[fld];
            }
        }
    }

#if DEBUG
    fclose(fp);
#endif
}