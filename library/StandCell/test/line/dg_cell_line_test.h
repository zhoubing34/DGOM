//
// Created by li12242 on 17/3/15.
//

#ifndef DGOM_DG_CELL_LINE_TEST_H
#define DGOM_DG_CELL_LINE_TEST_H

int dg_line_info_test(dg_cell *cell, int verbose);
int dg_line_node_test(dg_cell *cell, int verbose);
int dg_line_vand_matrix_test(dg_cell *cell, int verbose);
int dg_line_mass_matrix_test(dg_cell *cell, int verbose);
int dg_line_deri_matrix_test(dg_cell *cell, int verbose);
int dg_line_Fmask_test(dg_cell *cell, int verbose);
int dg_line_LIFT_test(dg_cell *quad, int verbose);
int dg_line_vert_proj_test(dg_cell *cell, int verbose);

#endif //DGOM_DG_CELL_LINE_TEST_H
