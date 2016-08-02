//
// Created by li12242 on 16/7/30.
//

#ifndef DGOM_SWEDRIVER2D_H
#define DGOM_SWEDRIVER2D_H

#include "PhysDomain/PhysDomain.h"
#include "LibUtilities/NetcdfLibrary.h"


typedef struct foo{
    double FinalTime;
    double dt;
    double gra;
    double dx;
    char   *casename;
    double **bot;
} SWESolver;


MultiReg2d*   SWEMesh2d(char *casename, SWESolver *solver, StdRegions2d *shape, int Mx, int My);
PhysDomain2d* SWEInit2d(char *casename, SWESolver *solver, MultiReg2d *mesh);
NcFile*       SWEOutput(PhysDomain2d *phys, SWESolver *solver);
void          StoreVar(NcFile *file, PhysDomain2d *phys, int outStep, double time);

#endif //DGOM_SWEDRIVER2D_H
