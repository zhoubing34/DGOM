#ifndef CONVECTION2D_H
#define CONVECTION2D_H

#include "PhysField/pf_phys.h"
#include "Utility/nc_library.h"

typedef enum{
    conv_rotational_convection = 1,
    conv_advection_diffusion = 2,
    conv_userset = 3,
} conv_caseid;

#define DIFF_THRESHOLD 1e-10

typedef struct {
    int caseid; ///< indicator of problem case
    char *casename; ///< case name
    sc_cellType celltype; ///< cell type
    int N; ///< order of polynomial
    int Ne; ///< number of elements in each coordinate
    double cfl; ///< CFL number
    double dt; ///< delta time
    double finaltime; ///< final time
    double u; ///< const flow rate on x direction
    double v; ///< const flow rate on y direction
    double viscosity; ///< diffusion parameter
    int isverbose; ///< parameters for debug
    nc_file *outfile; ///< output file
    double LDG_parameter[3]; ///< parameter C11,C12,C22 for LDG viscosity soler
} conv_solver2d;

#endif