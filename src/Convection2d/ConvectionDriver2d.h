#ifndef CONVECTION2D_H
#define CONVECTION2D_H

#include "PhysDomain/PhysDomain.h"
#include "NcOutput.h"

/* Mesh2d.c */
MultiReg2d* ReadMesh(StdRegions2d *shape, int Ne);
void PrintMesh ( MultiReg2d *mesh );

/* InitCondition.c */
double InitCondition(PhysDomain2d * phys, PhysDomain2d *flowRate);
void PrintPhys( PhysDomain2d *phys, char *name );

void ConvectionRun2d(PhysDomain2d *phys, PhysDomain2d *flowRate, Ncfile * outfile, double FinalTime, double dt);
void ConvectionRHS2d(PhysDomain2d *phys, PhysDomain2d *flowRate,
                     float frka, float frkb, float fdt);

/* Postprocess.c */
void Postprocess(PhysDomain2d *phys);

#endif