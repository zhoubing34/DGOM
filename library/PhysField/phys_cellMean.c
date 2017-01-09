//
// Created by li12242 on 17/1/9.
//

#include "phys_cellMean.h"

/**
 * @brief calculate the integral averaged cell value of each physical field
 * @param[in,out] phys physical field object
 */
void phys_cellMean(physField *phys){

    const int Nfield = phys->Nfield;
    const int K = phys->grid->K;
    const int Np = phys->cell->Np;

    register int k,n,fld,sk;

    real *f_Q = phys->f_Q;
    real *c_Q = phys->c_Q;
    real var[Np];
    for(k=0;k<K;k++){
        const double Area = phys->region->size[k];
        for(fld=0;fld<Nfield;fld++){
            sk = fld+(k*Np)*Nfield;
            for(n=0;n<Np;n++){
                var[n] = f_Q[sk]; sk+=Nfield;
            }
            c_Q[k*Nfield + fld] = mr_reg_integral(phys->region, k, var)/Area;
        }
    }

    return;
}