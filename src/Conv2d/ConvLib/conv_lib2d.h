//
// Created by li12242 on 17/4/9.
//

#ifndef DGOM_CONV_LIB2D_H
#define DGOM_CONV_LIB2D_H

#include "PhysField/dg_phys.h"
#include "PhysField/dg_phys_strong_surf_opt.h"
#include "PhysField/dg_phys_strong_vol_opt.h"
#include "PhysField/dg_phys_strong_LDG_opt.h"

typedef enum{
    CONV_HELP = 1,
    CONV_CREATE_INPUT = 2,
    CONV_RUN = 3,
} Conv_Run_Type;

typedef struct {
    dg_phys *phys;
    char filename[MAX_NAME_LENGTH]; ///< name of input file;
    double dt; ///< delta time;
    double outDt; ///< output time interval;
    double finaltime; ///< final time;
    NC_File *outfile; ///< output file;
    int vis_flag; ///< flag for adding viscosity;
} Conv_Solver;

/* conv_output */
void conv_setoutput();
void conv_putvar(dg_phys *phys, int timestep, double time);
/* conv_run */
void conv_run();
/* conv_rhs */
void conv_rhs(dg_phys *phys, dg_real frka, dg_real frkb, dg_real fdt);

#endif //DGOM_CONV_LIB2D_H
