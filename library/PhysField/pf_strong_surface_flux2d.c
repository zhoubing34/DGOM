//
// Created by li12242 on 12/27/16.
//

#include <StandCell/sc_stdcell.h>
#include "pf_strong_surface_flux2d.h"

#define DEBUG 1

/**
 * @brief calculate the surface integral of strong form
 * @details
 * \f[ M^{-1} \cdot M_{es} (E \cdot nx + G \cdot ny - Fhs) \f]
 * where the \f[ Fhs \f] is the numerical flux, \f[ F = (E, G) \f] is the flux term.
 *
 * @param[in,out] phys physField object
 * @param[in] slipwall_condition boundary condition for slip-wall
 * @param[in] nonslipwall_condition boundary condition for non-slip-wall
 * @param[in] nodal_flux function to calculate the flux term
 * @param[in] numerical_flux numerical flux function
 *
 * @note
 *
 */
void pf_strong_surface_flux2d
        (physField *phys,
         wall_condition_func slipwall_condition,  // slip wall condition
         wall_condition_func non_slipwall_condition, // non-slip wall condition
         nodal_flux_func nodal_flux, // flux term
         numerical_flux_func numerical_flux){ // numerical flux function

    const int K = phys->grid->K;
    const int Nfaces = phys->cell->Nfaces;
    const int Nfp = phys->cell->Nfp;
    const int Np = phys->cell->Np;
    const int Nfield = phys->Nfield;

    real *surfinfo = phys->surfinfo;
    real *f_Q = phys->f_Q;
    real *f_inQ = phys->f_inQ;
    real *f_ext = phys->f_ext;
    real *f_LIFT  = phys->cell->f_LIFT;

    real f_M[Nfield];
    real f_P[Nfield];
    real Eflux[Nfield], Gflux[Nfield], Fhs[Nfield];
    real flux_Q[Nfp*Nfaces*Nfield];

    register int k,m,n,fld,surfid=0;

    for(k=0; k<K; k++){
        for(m=0; m<Nfp*Nfaces; ++m){
            int  idM = (int)surfinfo[surfid++];
            int  idP = (int)surfinfo[surfid++];
            const real fsc = surfinfo[surfid++];
            const int  bstype = (int)surfinfo[surfid++];
            const real nx = surfinfo[surfid++];
            const real ny = surfinfo[surfid++];

            // local face values
            for(fld=0;fld<Nfield;fld++)
                f_M[fld] = f_Q[idM++];

            // adjacent nodal values
            switch (bstype){
                case INNERLOC:
                    for(fld=0;fld<Nfield;fld++)
                        f_P[fld] = f_Q[idP++];
                    break;
                case INNERBS:
                    for(fld=0;fld<Nfield;fld++)
                        f_P[fld] = f_inQ[idP++];
                    break;
                case SLIPWALL:
                    slipwall_condition(nx, ny, f_M, f_P); break;
                case NSLIPWALL:
                    non_slipwall_condition(nx, ny, f_M, f_P); break;
                default:
                    for(fld=0;fld<Nfield;fld++)
                        f_P[fld] = f_ext[idP++];
                    break;
            }

#if DEBUG
            nodal_flux(f_M, Eflux, Gflux);
            if( numerical_flux(nx, ny, f_M, f_P, Fhs) ){
                printf("pf_strong_surface_flux2d (%s): flux error at line %d\n"
                               "k=%d, m=%d\nbstype=%d,idP=%d\n",
                       __FILE__, __LINE__, k, m, bstype, idP);
                MPI_Abort(MPI_COMM_WORLD, -1);
            }
#else
            nodal_flux(f_M, Eflux, Gflux);
            numerical_flux(nx, ny, f_M, f_P, Fhs);
#endif
            real *t = flux_Q + m*Nfield;
            for(fld=0;fld<Nfield;fld++){
                t[fld] = fsc*( nx*Eflux[fld] + ny*Gflux[fld] - Fhs[fld] );
            }
        }

        for(n=0;n<Np;n++){
            const real *ptLIFT = f_LIFT + n*Nfp*Nfaces;
            real *f_rhsQ = phys->f_rhsQ + Nfield*(n+k*Np);
            for(m=0;m<Nfp*Nfaces;m++){
                const real L = ptLIFT[m];
                const real *t = flux_Q+m*Nfield;

                for(fld=0;fld<Nfield;fld++)
                    f_rhsQ[fld] += L*t[fld];
            }
        }
    }

    return;
}