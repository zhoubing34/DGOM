//
// Created by li12242 on 17/1/13.
//

#ifndef DGOM_SWE_MESH_H
#define DGOM_SWE_MESH_H

#include "swe_driver2d.h"

physField* swe_uniform_mesh(SWE_Solver *solver, int Mx, int My,
                            double xmin, double xmax, double ymin, double ymax);
physField* swe_file_mesh(SWE_Solver *solver, char *meshfile);
dg_real* swe_read_topography(SWE_Solver *solver, char *botfile);
dg_real* swe_flat_topography(SWE_Solver *solver);
dg_real* swe_parabolic_topography(SWE_Solver *solver);

#endif //DGOM_SWE_MESH_H
