#ifndef CONVECTION2D_H
#define CONVECTION2D_H

#include "PhysField/phys_physField.h"
#include "LibUtilities/nc_library.h"

typedef struct {
    int N; ///> order
    int Ne; ///> number of elements in each coordinate
    double cfl; ///> CFL number
    double dt; ///> delta time
    cellType type; ///> cell type
    double finalTime; ///> final time
    int isverbose; ///> parameters for debug
    nc_file *outfile; ///> output file
} conv_solver2d;

#endif