//
// Created by li12242 on 16/12/17.
//

#ifndef DGOM_SC_QUAD_TEST_H
#define DGOM_SC_QUAD_TEST_H

#include "StandCell/test/sc_test.h"

int sc_quadCoor_test(stdCell *quad, int verbose);
int sc_quadVandMatrix_test(stdCell *quad, int verbose);
int sc_quadMassMatrix_test(stdCell *quad, int verbose);
int sc_quadDeriMatrix_test(stdCell *quad, int verbose);
int sc_quadLIFT_test(stdCell *quad, int verbose);

#endif //DGOM_SC_QUAD_TEST_H
