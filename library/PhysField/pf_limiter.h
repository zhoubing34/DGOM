//
// Created by li12242 on 16/8/3.
//

#ifndef DGOM_SLOPELIMITER_H
#define DGOM_SLOPELIMITER_H

#include "PhysField/pf_phys.h"
void pf_slloc2d(physField *phys, double beta);
void pf_slop_limit(physField *phys, double beta);
void pf_vert_limit(physField *phys);

#endif //DGOM_SLOPELIMITER_H
