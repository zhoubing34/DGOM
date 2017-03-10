//
// Created by li12242 on 12/19/16.
//

#ifndef DGOM_MR_REG_TEST_H
#define DGOM_MR_REG_TEST_H

#include "multiregion_test.h"
int dg_region_node_test(dg_region *reg, int verbose);
int dg_region_volume_factor_test(dg_region *reg, int verbose);
int dg_region_face_factor_test(dg_region *reg, int verbose);
int dg_region_scale_test(dg_region *reg, int verbose);

#endif //DGOM_MR_REG_TEST_H
