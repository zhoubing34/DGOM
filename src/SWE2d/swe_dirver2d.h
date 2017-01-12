//
// Created by li12242 on 16/7/30.
//

#ifndef DGOM_SWEDRIVER2D_H
#define DGOM_SWEDRIVER2D_H

#include "PhysField/phys_physField.h"
#include "LibUtilities/nc_library.h"
#include "LibUtilities/LibUtilities.h"

typedef enum {
    swe_dambreakwet = 1,
    swe_dambreakdry = 2,
    swe_parabolicbowl = 3,
    swe_userset = 4,
}casename;

typedef struct{
    int caseid; ///< case id
    /* standard cell */
    sc_cellType celltype; ///< type id of standard cell
    int N; ///< degree of polynomial
    int Mx; ///< No. of cells on x direction
    int My; ///< No. of cells on y direction
    char *meshname; ///< mesh name
    /* time */
    double cfl; ///< cfl number
    double dt; ///< minimum time step
    double ftime; ///< final time
    /* physical parameters */
    double gra; ///< gravity acceleration
    char *botname; ///< file name of bottom topography
    double *bot; ///< bottom topography
    double hcrit; ///< minimum water depth
    double roughness; ///< manning roughness coefficient for friction term
    /* output */
    nc_file *outfile; ///< netCDF output file
    /* LDG parameter */
    double LDG_parameter[3];
}swe_solver;


#endif //DGOM_SWEDRIVER2D_H
