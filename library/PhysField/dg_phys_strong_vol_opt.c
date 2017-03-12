//
// Created by li12242 on 12/26/16.
//

#include "dg_phys_strong_vol_opt.h"

#define DEBUG 0

/**
 * @brief RHS of volume flux term F = (E, G) integral for 2-d problem
 * @details
 *
 * @param[in,out] phys
 * @param[in] nodal_flux
 *
 * @note
 * The field of variable `phys->rhs` will be initialized and assigned to new values.
 */
void dg_phys_strong_vol_opt2d(dg_phys *phys, Nodal_Flux_Fun nodal_flux){

    const int K = dg_grid_K(phys->grid);
    const int Np = dg_cell_Np(phys->cell);
    const int Nfield = dg_phys_Nfield(phys);

    dg_real *f_Q = phys->f_Q;
    dg_real *f_rhsQ = phys->f_rhsQ;
    dg_real *f_Dr = phys->cell->f_Dr;
    dg_real *f_Ds = phys->cell->f_Ds;

    register unsigned int k,n,m,fld,rhsid=0;
    dg_real **drdx_p = phys->region->drdx;
    dg_real **drdy_p = phys->region->drdy;
    dg_real **dsdx_p = phys->region->dsdx;
    dg_real **dsdy_p = phys->region->dsdy;

    dg_real Eflux[Np*Nfield], Gflux[Np*Nfield], rhs[Nfield];

    for(k=0;k<K;k++){
        dg_real *var = f_Q + k*Np*Nfield; // variable in k-th element

        // calculate flux term
        for(n=0;n<Np;n++){
            // calculate flux term on the n-th point
            nodal_flux(var+n*Nfield, Eflux+n*Nfield, Gflux+n*Nfield);
#if DEBUG
            if(!phys->grid->procid){
                for(fld = 0; fld<Nfield; fld ++)
                    printf("k=%d, n=%d, fld=%d, Eflux=%f, Gflux=%f\n", k, n, fld, Eflux[n*Nfield+fld], Gflux[n*Nfield+fld]);
            }
#endif
        }

        for(n=0;n<Np;++n){ // rhs for n-th point

            const dg_real *ptDr = f_Dr+n*Np; // n-th row of Dr
            const dg_real *ptDs = f_Ds+n*Np; // n-th row of Ds

            const dg_real drdx = drdx_p[k][n]; // volume geometry for n-th point
            const dg_real drdy = drdy_p[k][n]; // volume geometry for n-th point
            const dg_real dsdx = dsdx_p[k][n]; // volume geometry for n-th point
            const dg_real dsdy = dsdy_p[k][n]; // volume geometry for n-th point

            // initialize rhs
            for(m=0;m<Nfield;m++){
                rhs[m] = 0;
            }
#if DEBUG
            if(!phys->grid->procid)
                printf("k=%d, n=%d, drdx=%f, drdy=%f, dsdx=%f, dsdy=%f\n", k, n, drdx, drdy, dsdx, dsdy);
#endif
            for(m=0;m<Np;++m){
                const dg_real dr = ptDr[m]; // m-th column for Dr
                const dg_real ds = ptDs[m]; // m-th column for Ds
                const dg_real dx = drdx*dr+dsdx*ds;
                const dg_real dy = drdy*dr+dsdy*ds;

                const dg_real *eflux = Eflux + m*Nfield;
                const dg_real *gflux = Gflux + m*Nfield;

#if DEBUG
                if(!phys->grid->procid)
                    printf("k=%d, m=%d, Dr=%f, Ds=%f, dx=%f, dy=%f\n", k, m, dr, ds, dx, dy);
#endif

                for(fld=0;fld<Nfield;fld++){ // m-th point over all field
                    rhs[fld] += -(dx*eflux[fld] + dy*gflux[fld]);
                }
            }

            for(fld=0;fld<Nfield;fld++){
                f_rhsQ[rhsid++] = rhs[fld];
            }
        }
    }
}