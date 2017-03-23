//
// Created by li12242 on 17/2/20.
//

#ifndef DGOM_PF_OPENBC_H
#define DGOM_PF_OPENBC_H

#include "Utility/nc_library.h"
#include "Utility/utility.h"
#include "dg_phys_info.h"

typedef enum{
    INTERP_LINEAR = 0,
} dg_phys_obc_interp_type;

typedef struct dg_phys_obc{
    dg_phys_info *info;
    NC_File *file;
    dg_real *f_extQ;

    int Ntime; ///< number of time steps in obc file;
    int Nvert; ///< number of vertex in obc file;
    double *time; ///< time vector in obc file;
    int *vert; ///< index of vertex in obc file;
    dg_phys_obc_interp_type interp_type; ///< interpolation type;

    void (*add_obc)(struct dg_phys_obc *obc, char *filename);
    void (*update_obc)(struct dg_phys_obc *obc, double time);
}dg_phys_obc;

dg_phys_obc* dg_phys_obc_create(dg_phys_info *phys_info);
void dg_phys_obc_free(dg_phys_obc *phys_obc);

#endif //DGOM_PF_OPENBC_H
