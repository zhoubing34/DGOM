//
// Created by li12242 on 16/12/30.
//

#include "pf_add_LDG_solver.h"
#include "pf_phys.h"

void pf_add_LDG_solver(physField *phys){

    phys->viscosity = (pf_LDG_solver*) calloc(1, sizeof(pf_LDG_solver));

    const int Nfield = phys->Nfield;
    const int K = phys->grid->K;
    const int Np = phys->cell->Np;

    phys->viscosity->px_Q = (dg_real *) calloc((size_t) Nfield*K*Np, sizeof(dg_real));
    phys->viscosity->py_Q = (dg_real *) calloc((size_t) Nfield*K*Np, sizeof(dg_real));
    phys->viscosity->vis_Q = (dg_real *) calloc((size_t) Nfield*K*Np, sizeof(dg_real));

    const int parallNodeNum = phys->parallNodeNum;

    phys->viscosity->px_inQ = (dg_real *) calloc((size_t) parallNodeNum, sizeof(dg_real));
    phys->viscosity->py_inQ = (dg_real *) calloc((size_t) parallNodeNum, sizeof(dg_real));
    phys->viscosity->px_outQ = (dg_real *) calloc((size_t) parallNodeNum, sizeof(dg_real));
    phys->viscosity->py_outQ = (dg_real *) calloc((size_t) parallNodeNum, sizeof(dg_real));

    return;
}

void pf_delete_LDG_solver(physField *phys){

    free(phys->viscosity->px_Q);
    free(phys->viscosity->py_Q);
    free(phys->viscosity->vis_Q);

    free(phys->viscosity->px_inQ);
    free(phys->viscosity->py_inQ);

    free(phys->viscosity->px_outQ);
    free(phys->viscosity->py_outQ);

    free(phys->viscosity);
}