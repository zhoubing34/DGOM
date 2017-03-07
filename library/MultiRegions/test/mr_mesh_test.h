//
// Created by li12242 on 12/21/16.
//

#ifndef DGOM_MR_MESH_TEST_H
#define DGOM_MR_MESH_TEST_H

#include "mr_test.h"
int mr_mesh_connet_test(dg_mesh *mesh, int verbose);
int mr_mesh_parallel_test(dg_mesh *mesh, int verbose);
int mr_mesh_bc_test(dg_mesh *mesh, int verbose);

#endif //DGOM_MR_MESH_TEST_H
