//
// Created by li12242 on 17/2/20.
//

#ifndef DGOM_PF_OPENBC_H
#define DGOM_PF_OPENBC_H

#include "dg_phys.h"

typedef enum{
    time_interp_linear=0
} time_interp_method;

void pf_set_openbc(physField *phys, double timeloc, time_interp_method method);

#endif //DGOM_PF_OPENBC_H
