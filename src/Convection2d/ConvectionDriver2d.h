#ifndef CONVECTION2D_H
#define CONVECTION2D_H

#include "PhysDomain/PhysDomain.h"
#include "LibUtilities/NetcdfLibrary.h"

/* Mesh2d.c */
MultiReg2d* ReadMesh(StdRegions2d *shape, int Ne);

/* InitCondition.c */
double InitCondition(PhysDomain2d * phys, PhysDomain2d *flowRate);
void PrintPhys( PhysDomain2d *phys, char *name );

void ConvectionRun2d(PhysDomain2d *phys, PhysDomain2d *flowRate, NcFile * outfile, double FinalTime, double dt);
void ConvectionRHS2d(PhysDomain2d *phys, PhysDomain2d *flowRate,
                     float frka, float frkb, float fdt);

/* Initialize output */
NcFile* SetupOutput(MultiReg2d *mesh, char* casename);
void PutVar(NcFile *file, int outStep, double time, PhysDomain2d* phys);

/* Postprocess.c */
void Postprocess(PhysDomain2d *phys);

#endif