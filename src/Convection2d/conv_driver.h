#ifndef CONVECTION2D_H
#define CONVECTION2D_H

#include "Utility/nc_library.h"
#include "Utility/arg_section.h"
#include "MultiRegions/Grid/dg_grid_reader.h"
#include "MultiRegions/Mesh/dg_mesh_fetch_buffer.h"
#include "PhysField/dg_phys.h"
#include "PhysField/dg_phys_strong_surf_opt.h"
#include "PhysField/dg_phys_strong_vol_opt.h"

#define SEC_NUM 6
#define HEADEND   "[============]"
#define HEADLINE  "[------------]"

typedef enum{
    conv_help = 1,
    conv_preprocss = 2,
    conv_execute = 3,
} Conv_Run_Type;

typedef enum{
    conv_rotational_convection = 1,
    conv_advection_diffusion = 2,
    conv_userset = 3,
} Conv_Case_Type;

#define INPUT_FILE_NAME "conv2d_paramter.inc"

typedef struct {
    dg_phys *phys;
    double dt; ///< delta time
    double finaltime; ///< final time
    nc_file *outfile; ///< output file
} Conv_Solver;

/* conv_input */
Conv_Case_Type conv_input(int argc, char *agrv[]);
void conv_arg_section_free(arg_section **section_p);
arg_section** conv_read_input();
/* conv_init */
void conv_init(Conv_Case_Type case_type);
/* conv_output */
void conv_setoutput();
void conv_putvar(dg_phys *phys, int timestep, double time);
/* conv_run */
void conv_run();
/* conv_rhs */
void conv_rhs(dg_phys *phys, dg_real frka, dg_real frkb, dg_real fdt);
/* conv_extsol */
void conv_normerr(Conv_Case_Type case_type);


#endif