//
// Created by li12242 on 16/12/16.
//

#ifndef DGOM_SC_TRI_TEST_H
#define DGOM_SC_TRI_TEST_H

#include "StandCell/test/sc_test_main.h"

int sc_triCoor_test(dg_cell *tri, int verbose);
int sc_triVandMatrix_test(dg_cell *tri, int verbose);
int sc_triMassMatrix_test(dg_cell *tri, int verbose);
int sc_triDeriMatrix_test(dg_cell *tri, int verbose);
int sc_triLIFT_test(dg_cell *tri, int verbose);
int sc_triVertProj_test(dg_cell *tri, int verbose);

#endif //DGOM_SC_TRI_TEST_H
