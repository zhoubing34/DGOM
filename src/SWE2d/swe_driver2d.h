//
// Created by li12242 on 16/7/30.
//

#ifndef DGOM_SWEDRIVER2D_H
#define DGOM_SWEDRIVER2D_H

#include "PhysField/dg_phys.h"
#include "Utility/nc_library.h"
#include "Utility/utility.h"

typedef enum {
    swe_dambreakwet = 1,
    swe_dambreakdry = 2,
    swe_parabolicbowl = 3,
    swe_userset = 4,
}SWE_Run_Type;

typedef struct{
    /* physical field */
    dg_phys *phys; ///< physical field
    double ftime; ///< final time
    /* physical parameters */
    double gra; ///< gravity acceleration
    double hcrit; ///< minimum water depth
    double *bot; ///< bottom topography
    double *rough; ///< manning roughness coefficient for friction term
    /* output file */
    NC_File *ncfile; ///< netCDF output file
}SWE_Solver;


#endif //DGOM_SWEDRIVER2D_H
