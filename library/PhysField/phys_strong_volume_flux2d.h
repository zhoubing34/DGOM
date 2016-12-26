//
// Created by li12242 on 12/26/16.
//

#ifndef DGOM_PHYS_STRONG_VOLUME_FLUX2D_H
#define DGOM_PHYS_STRONG_VOLUME_FLUX2D_H

#include "phs_physField.h"

// nodal flux term function for Nfield variables
typedef int( *nodal_flux_func)(real *var, real *Eflux, real *Gflux);

void phys_strong_volume_flux2d(physField *phys, nodal_flux_func nodal_flux);

#endif //DGOM_PHYS_STRONG_VOLUME_FLUX2D_H
