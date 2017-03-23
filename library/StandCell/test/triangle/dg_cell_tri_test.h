//
// Created by li12242 on 16/12/16.
//

#ifndef DGOM_SC_TRI_TEST_H
#define DGOM_SC_TRI_TEST_H

#include "StandCell/test/dg_cell_test_main.h"

int dg_tri_info_test(dg_cell *cell, int verbose);
int dg_tri_nood_test(dg_cell *tri, int verbose);
int dg_tri_vand_matrix_test(dg_cell *tri, int verbose);
int dg_tri_mass_matrix_test(dg_cell *tri, int verbose);
int dg_tri_deri_matrix_test(dg_cell *tri, int verbose);
int dg_tri_Fmask_test(dg_cell *tri, int verbose);
int dg_tri_LIFT_test(dg_cell *tri, int verbose);
int dg_tri_vert_proj_test(dg_cell *tri, int verbose);

#endif //DGOM_SC_TRI_TEST_H
