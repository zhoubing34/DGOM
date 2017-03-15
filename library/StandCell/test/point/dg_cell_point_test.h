//
// Created by li12242 on 17/3/15.
//

#ifndef DGOM_DG_CELL_POINT_TEST_H
#define DGOM_DG_CELL_POINT_TEST_H

int dg_point_info_test(dg_cell *cell, int verbose);
int dg_point_node_test(dg_cell *cell, int verbose);
int dg_point_vand_matrix_test(dg_cell *cell, int verbose);
int dg_point_mass_matrix_test(dg_cell *cell, int verbose);
int dg_point_deri_matrix_test(dg_cell *cell, int verbose);
int dg_point_Fmask_test(dg_cell *cell, int verbose);
int dg_point_LIFT_test(dg_cell *quad, int verbose);

#endif //DGOM_DG_CELL_POINT_TEST_H
