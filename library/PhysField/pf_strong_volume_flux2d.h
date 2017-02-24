//
// Created by li12242 on 12/26/16.
//

#ifndef DGOM_PHYS_STRONG_VOLUME_FLUX2D_H
#define DGOM_PHYS_STRONG_VOLUME_FLUX2D_H

#include "pf_phys.h"

// nodal flux term function for Nfield variables
typedef int( *nodal_flux_func)(dg_real *var, dg_real *Eflux, dg_real *Gflux);

void pf_strong_volume_flux2d(physField *phys, nodal_flux_func nodal_flux);

#endif //DGOM_PHYS_STRONG_VOLUME_FLUX2D_H
