//
// Created by li12242 on 16/12/30.
//

#include "phys_add_LDG_solver.h"

void phys_add_LDG_solver(physField *phys){

    const int Nfield = phys->Nfield;
    const int K = phys->grid->K;
    const int Np = phys->cell->Np;

    phys->viscosity.px_Q = (real *) calloc((size_t) Nfield*K*Np, sizeof(real));
    phys->viscosity.py_Q = (real *) calloc((size_t) Nfield*K*Np, sizeof(real));

    const int parallNodeNum = phys->parallNodeNum;

    phys->viscosity.px_inQ = (real *) calloc((size_t) parallNodeNum, sizeof(real));
    phys->viscosity.py_inQ = (real *) calloc((size_t) parallNodeNum, sizeof(real));
    phys->viscosity.px_outQ = (real *) calloc((size_t) parallNodeNum, sizeof(real));
    phys->viscosity.py_outQ = (real *) calloc((size_t) parallNodeNum, sizeof(real));

    phys->viscosity.vis_Q = (real *) calloc((size_t) Nfield*K*Np, sizeof(real));

    return;
}

void phys_delete_LDG_solver(physField *phys){

    free(phys->viscosity.px_Q);
    free(phys->viscosity.py_Q);

    free(phys->viscosity.px_inQ);
    free(phys->viscosity.py_inQ);

    free(phys->viscosity.px_outQ);
    free(phys->viscosity.px_outQ);

    free(phys->viscosity.vis_Q);
}