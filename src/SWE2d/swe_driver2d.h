//
// Created by li12242 on 16/7/30.
//

#ifndef DGOM_SWEDRIVER2D_H
#define DGOM_SWEDRIVER2D_H

#include "PhysField/pf_phys.h"
#include "Utility/nc_library.h"
#include "Utility/LibUtilities.h"

typedef enum {
    swe_dambreakwet = 1,
    swe_dambreakdry = 2,
    swe_parabolicbowl = 3,
    swe_userset = 4,
}swe_case;

typedef struct{
    swe_case caseid; ///< case id
    char *casename;
    /* standard cell */
    sc_cellType celltype; ///< type id of standard cell
    int N; ///< degree of polynomial
    /* physical field */
    physField *phys;
    /* time */
    double cfl; ///< cfl number
    //double dt; ///< minimum time step
    double ftime; ///< final time
    /* physical parameters */
    double gra; ///< gravity acceleration
    double *bot; ///< bottom topography
    double hcrit; ///< minimum water depth
    double roughness; ///< manning roughness coefficient for friction term
    /* output file */
    nc_file *outfile; ///< netCDF output file
    /* LDG parameter */
    double c12;
}swe_solver;


#endif //DGOM_SWEDRIVER2D_H
