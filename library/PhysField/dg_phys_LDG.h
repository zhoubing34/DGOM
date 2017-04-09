//
// Created by li12242 on 16/12/30.
//

#ifndef DGOM_PHYS_ADD_VISCOSITY_SOLVER_H
#define DGOM_PHYS_ADD_VISCOSITY_SOLVER_H

#include "dg_phys_info.h"

typedef struct dg_phys_LDG{
    dg_phys_info *info; ///< pointer to dg_phys_info structure;
    dg_real *px, *py, *pz; ///< auxiliary variale in LDG methods;
    dg_real *px_recv, *py_recv, *pz_recv; ///< buffers for receiving data;
    dg_real *sqrt_miux, *sqrt_miuy, *sqrt_miuz; ///< sqrt root of viscosity;
} dg_phys_LDG;

dg_phys_LDG* dg_phys_LDG_create(dg_phys_info *info);
void dg_phys_LDG_free(dg_phys_LDG *ldg);

#define dg_phys_ldg_px(ldg) (ldg->px)
#define dg_phys_ldg_py(ldg) (ldg->py)
#define dg_phys_ldg_pz(ldg) (ldg->pz)
#define dg_phys_ldg_sqrt_miux(ldg) (ldg->sqrt_miux)
#define dg_phys_ldg_sqrt_miuy(ldg) (ldg->sqrt_miuy)
#define dg_phys_ldg_sqrt_miuz(ldg) (ldg->sqrt_miuz)
#define dg_phys_ldg_x_recv(ldg) (ldg->px_recv)
#define dg_phys_ldg_y_recv(ldg) (ldg->py_recv)
#define dg_phys_ldg_z_recv(ldg) (ldg->pz_recv)

#endif //DGOM_PHYS_ADD_VISCOSITY_SOLVER_H
