//
// Created by li12242 on 12/27/16.
//

#include <StandCell/sc_stdcell.h>
#include "phys_strong_surface_flux2d.h"

#define hand_err(ret) do{\
if( ret ){\
    printf("PhysField (%s): flux error at line %d\n", __FILE__, __LINE__);\
    MPI_Abort(MPI_COMM_WORLD, -1);\
}\
}while(0) \

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
void phys_strong_surface_integral2d
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

    register int k,m,n,fld, surfid=0;

    real *surfinfo = phys->surfinfo;
    real *f_Q = phys->f_Q;
    real *f_inQ = phys->f_inQ;
    real *f_ext = phys->f_ext;
    real *f_LIFT  = phys->cell->f_LIFT;

    real varM[Nfield], varP[Nfield];
    real Eflux[Nfield], Gflux[Nfield], Fhs[Nfield];
    real dflux[Nfp*Nfaces*Nfield];

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
                varM[fld] = f_Q[idM++];

            // adjacent nodal values
            switch (bstype){
                case INNERLOC:
                    for(fld=0;fld<Nfield;fld++)
                        varP[fld] = f_Q[idP++];
                    break;
                case INNERBS:
                    for(fld=0;fld<Nfield;fld++)
                        varP[fld] = f_inQ[idP++];
                    break;
                case SLIPWALL:
                    slipwall_condition(nx, ny, varM, varP); break;
                case NSLIPWALL:
                    non_slipwall_condition(nx, ny, varM, varP); break;
                default:
                    for(fld=0;fld<Nfield;fld++)
                        varP[fld] = f_ext[idP++];
                    break;
            }

            nodal_flux(varM, Eflux, Gflux);
            numerical_flux(nx, ny, varM, varP, Fhs);
//            hand_err( nodal_flux(varM, Eflux, Gflux) );
//            hand_err( numerical_flux(nx, ny, varM, varP, Fhs) );

            real *flux_Q = dflux + m*Nfield;
            for(fld=0;fld<Nfield;fld++){
                flux_Q[fld] = fsc*( nx*Eflux[fld] + ny*Gflux[fld] - Fhs[fld] );
            }
        }

        for(n=0;n<Np;n++){
            const real *ptLIFT = f_LIFT + n*Nfp*Nfaces;
            for(m=0;m<Nfp*Nfaces;m++){
                const real L = ptLIFT[m];
                real *flux_Q = dflux+m*Nfield;
                real *f_rhsQ = phys->f_rhsQ + Nfield*(n+k*Np);

                for(fld=0;fld<Nfield;fld++)
                    f_rhsQ[fld] += L*flux_Q[fld];
            }
        }
    }

    return;
}