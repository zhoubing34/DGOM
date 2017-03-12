//
// Created by li12242 on 12/26/16.
//

#ifndef DGOM_PHYS_STRONG_VOLUME_FLUX2D_H
#define DGOM_PHYS_STRONG_VOLUME_FLUX2D_H

#include "dg_phys.h"

// nodal flux term function for Nfield variables
typedef int( *Nodal_Flux_Fun)(dg_real *var, dg_real *Eflux, dg_real *Gflux);
void dg_phys_strong_vol_opt2d(dg_phys *phys, Nodal_Flux_Fun nodal_flux);

#endif //DGOM_PHYS_STRONG_VOLUME_FLUX2D_H
