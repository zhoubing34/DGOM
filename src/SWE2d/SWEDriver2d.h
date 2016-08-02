//
// Created by li12242 on 16/7/30.
//

#ifndef DGOM_SWEDRIVER2D_H
#define DGOM_SWEDRIVER2D_H

#include "PhysDomain/PhysDomain.h"
#include "LibUtilities/NetcdfLibrary.h"


typedef struct foo{
    double gra;         /* gravity acceleration */
    double dx;          /* minimum grid length */
    double **bot;       /* bottom topography */
    double FinalTime;   /* final time */
    double hcrit;       /* minimum water depth */
    char   *casename;   /* casename */
} SWESolver;


MultiReg2d*   SWEMesh2d(char *casename, SWESolver *solver, StdRegions2d *shape, int Mx, int My);
PhysDomain2d* SWEInit2d(char *casename, SWESolver *solver, MultiReg2d *mesh);
NcFile*       SWEOutput(PhysDomain2d *phys, SWESolver *solver);
void          StoreVar(NcFile *file, PhysDomain2d *phys, int outStep, double time);

double        SWERun2d(PhysDomain2d *phys, SWESolver *solver, NcFile *outfile);
void          SWERHS2d(PhysDomain2d *phys, SWESolver *solver,
                       const real fa, const real fb, const real ft);

double        SWEPredictDt(PhysDomain2d *phys, SWESolver *solver, double CFL);
void          SWEFlux2d(PhysDomain2d *phys, SWESolver *solver,
                        real *Qk, real *Eflx, real *Gflx);
void          SWESource(PhysDomain2d *phys, SWESolver *solver,
                        int k, real *vgeo, real *Qk, real *Sour);
void          SWENumFlux2d(SWESolver *solver, real nx, real ny,
                           real hM, real hP, real qxM, real qxP, real qyM, real qyP,
                           real *Fhs, real *Fqxs, real *Fqys);

#endif //DGOM_SWEDRIVER2D_H
