#ifndef SWE_OUTPUT_H
#define SWE_OUTPUT_H

#include "swe_driver2d.h"

nc_file* swe_output(swe_solver *solver);
void swe_save_var(swe_solver *solver, int tstep, double time);

#endif