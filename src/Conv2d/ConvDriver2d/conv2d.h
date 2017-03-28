#ifndef CONVECTION2D_H
#define CONVECTION2D_H

#include "Utility/nc_library.h"
#include "Utility/arg_section.h"
#include "MultiRegions/Grid/dg_grid_reader.h"
#include "PhysField/dg_phys.h"
#include "PhysField/dg_phys_strong_surf_opt.h"
#include "PhysField/dg_phys_strong_vol_opt.h"

#define SEC_NUM 6
#define HEADEND   "[============]"
#define HEADLINE  "[------------]"

typedef enum{
    conv_help = 1,
    conv_create_input = 2,
    conv_execute = 3,
} Conv_Run_Type;

/** case ID for conv2d_st problem */
typedef enum{
    conv_rotational_convection = 1,
    conv_advection_diffusion = 2,
} Conv_St_Case;

typedef struct {
    dg_phys *phys;
    char filename[MAX_NAME_LENGTH]; ///< name of input file;
    double dt; ///< delta time;
    double finaltime; ///< final time;
    NC_File *outfile; ///< output file;
} Conv_Solver;

/* conv_input */
void conv_input(int argc, char *agrv[]);
void conv_arg_section_free(arg_section **section_p);
arg_section** conv_read_inputfile(char *);
/* conv_init */
void conv_init();
/* conv_output */
void conv_setoutput();
void conv_putvar(dg_phys *phys, int timestep, double time);
/* conv_run */
void conv_run();
/* conv_rhs */
void conv_rhs(dg_phys *phys, dg_real frka, dg_real frkb, dg_real fdt);
#endif