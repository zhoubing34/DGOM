#ifndef SWE_OUTPUT_H
#define SWE_OUTPUT_H

#include "swe_driver2d.h"

nc_file* swe_output(SWE_Solver *solver);
void swe_save_var(SWE_Solver *solver, int tstep, double time);

#endif