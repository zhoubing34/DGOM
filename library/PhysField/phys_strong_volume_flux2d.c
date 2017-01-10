//
// Created by li12242 on 12/26/16.
//

#include "phys_strong_volume_flux2d.h"

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
void phys_strong_volume_flux2d(physField *phys, nodal_flux_func nodal_flux){
    const int K = phys->grid->K;
    const int Np = phys->cell->Np;
    const int Nfield = phys->Nfield;

    real *f_Q = phys->f_Q;
    real *f_rhsQ = phys->f_rhsQ;
    real *f_Dr = phys->cell->f_Dr;
    real *f_Ds = phys->cell->f_Ds;
    real *vgeo = phys->vgeo;

    register unsigned int k,n,m,fld,geoid=0, rhsid=0;

    real Eflux[Np*Nfield], Gflux[Np*Nfield], rhs[Nfield];

    for(k=0;k<K;k++){
        real *var = f_Q + k*Np*Nfield; // variable in k-th element

        // calculate flux term
        for(n=0;n<Np;n++){
            // calculate flux term on the n-th point
            nodal_flux(var + n*Nfield, Eflux + n*Nfield, Gflux + n*Nfield);
//            if( nodal_flux(var + n*Nfield, Eflux + n*Nfield, Gflux + n*Nfield) ){
//                fprintf(stderr ,"PhysField (%s): flux error at k=%d, n=%d\n", __FILE__, k, n);
//                MPI_Abort(MPI_COMM_WORLD, -1);
//            }
#if DEBUG
            if(!phys->grid->procid){
                for(fld = 0; fld<Nfield; fld ++)
                    printf("k=%d, n=%d, fld=%d, Eflux=%f, Gflux=%f\n", k, n, fld, Eflux[n*Nfield+fld], Gflux[n*Nfield+fld]);
            }
#endif
        }

        for(n=0;n<Np;++n){ // rhs for n-th point

            const real *ptDr = f_Dr+n*Np; // n-th row of Dr
            const real *ptDs = f_Ds+n*Np; // n-th row of Ds

            const real drdx = vgeo[geoid++]; // volume geometry for n-th point
            const real drdy = vgeo[geoid++]; // volume geometry for n-th point
            const real dsdx = vgeo[geoid++]; // volume geometry for n-th point
            const real dsdy = vgeo[geoid++]; // volume geometry for n-th point

            // initialize rhs
            for(m=0;m<Nfield;m++){
                rhs[m] = 0;
            }
#if DEBUG
            if(!phys->grid->procid)
                printf("k=%d, n=%d, drdx=%f, drdy=%f, dsdx=%f, dsdy=%f\n", k, n, drdx, drdy, dsdx, dsdy);
#endif
            for(m=0;m<Np;++m){
                const real dr = ptDr[m]; // m-th column for Dr
                const real ds = ptDs[m]; // m-th column for Ds
                const real dx = drdx*dr+dsdx*ds;
                const real dy = drdy*dr+dsdy*ds;

                const real *eflux = Eflux + m*Nfield;
                const real *gflux = Gflux + m*Nfield;

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