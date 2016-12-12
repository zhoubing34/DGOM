//
// Created by li12242 on 16/7/30.
//

#ifndef DGOM_SWEDRIVER2D_H
#define DGOM_SWEDRIVER2D_H

#include "PhysDomain/PhysDomain.h"
#include "LibUtilities/NetcdfLibrary.h"
#include "PhysDomain/SlopeLimiter.h"
#include "LibUtilities/LibUtilities.h"

typedef struct{
    double gra;         /* gravity acceleration */
    double dx;          /* minimum grid length */
    double **bot;       /* bottom topography */
    double FinalTime;   /* final time */
    double hcrit;       /* minimum water depth */
    double dtmin;       /* minimum time step */
    char   *casename;   /* casename */
} SWE_Solver2d;


MultiReg2d*   SWE_Mesh2d(char **argv, SWE_Solver2d *solver);
PhysDomain2d* SWE_Init2d(char **argv, SWE_Solver2d *solver, MultiReg2d *mesh);
NcFile*       SWE_SetNcOutput2d(PhysDomain2d *phys, SWE_Solver2d *solver);
void          SWE_StoreVar2d(NcFile *file, PhysDomain2d *phys, int outStep, double time);

double        SWE_Run2d(PhysDomain2d *phys, SWE_Solver2d *solver, NcFile *outfile);
void          SWE_RHS2d(PhysDomain2d *phys, SWE_Solver2d *solver,
                        const real fa, const real fb, const real ft);

double        SWE_PredictDt2d(PhysDomain2d *phys, SWE_Solver2d *solver, double CFL);
void          SWE_ElementalFlux2d(PhysDomain2d *phys, SWE_Solver2d *solver,
                                  real *Qk, real *Eflx, real *Gflx);
void          SWE_NodalFlux2d(SWE_Solver2d *solver,
                              real h, real qx, real qy,
                              real *Eh, real *Eqx, real *Eqy,
                              real *Gh, real *Gqx, real *Gqy);
void          SWE_ElementalSource2d(PhysDomain2d *phys, SWE_Solver2d *solver,
                                    int k, real *vgeo, real *Qk, real *soureTerm);
void          SWE_NodalNumFlux2d(SWE_Solver2d *solver, real nx, real ny,
                                 real hM, real hP, real qxM, real qxP, real qyM, real qyP,
                                 real *Fhs, real *Fqxs, real *Fqys);
void          SWE_PositivePreserving2d(PhysDomain2d *phys, SWE_Solver2d *solver);

#endif //DGOM_SWEDRIVER2D_H
