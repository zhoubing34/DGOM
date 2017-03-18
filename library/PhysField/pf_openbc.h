//
// Created by li12242 on 17/2/20.
//

#ifndef DGOM_PF_OPENBC_H
#define DGOM_PF_OPENBC_H

#include "dg_phys.h"
#include "Utility/nc_library.h"
typedef enum{
    INTERP_LINEAR=0,
} dg_phys_obc_interp_type;

typedef struct dg_phys_obc{
    nc_file *file;
    dg_real *f_extQ;
    dg_phys_obc_interp_type interp_type;
}dg_phys_obc;

void pf_set_openbc(physField *phys, double timeloc, time_interp_method method);

#endif //DGOM_PF_OPENBC_H
