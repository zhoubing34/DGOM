#ifndef DGOM_CONV_OUTPUT
#define DGOM_CONV_OUTPUT

#include "PhysField/phys_physField.h"
void conv_setoutput(physField *phys);
void conv_putvar(physField *phys, int timestep, double time);

#endif