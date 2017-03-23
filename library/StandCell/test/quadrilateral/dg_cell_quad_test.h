//
// Created by li12242 on 16/12/17.
//

#ifndef DGOM_SC_QUAD_TEST_H
#define DGOM_SC_QUAD_TEST_H

#include "StandCell/test/dg_cell_test_main.h"

int dg_quad_info_test(dg_cell *cell, int verbose);
int dg_quad_nood_test(dg_cell *quad, int verbose);
int dg_quad_vand_matrix_test(dg_cell *quad, int verbose);
int dg_quad_mass_matrix_test(dg_cell *quad, int verbose);
int dg_quad_deri_matrix_test(dg_cell *quad, int verbose);
int dg_quad_Fmask_test(dg_cell *cell, int verbose);
int dg_quad_LIFT_test(dg_cell *quad, int verbose);
int dg_quad_vert_proj_test(dg_cell *cell, int verbose);

#endif //DGOM_SC_QUAD_TEST_H
