//
// Created by li12242 on 16/7/30.
//

#ifndef DGOM_SWEDRIVER2D_H
#define DGOM_SWEDRIVER2D_H

#include "PhysDomain/PhysDomain.h"


typedef struct foo{
    double FinalTime;
    double dt;
    double dx;
    char   *casename;
} Solver;


MultiReg2d* SWEMesh2d(char *casename, StdRegions2d *shape, int Mx, int My);

#endif //DGOM_SWEDRIVER2D_H
