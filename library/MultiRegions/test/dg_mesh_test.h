//
// Created by li12242 on 12/21/16.
//

#ifndef DGOM_MR_MESH_TEST_H
#define DGOM_MR_MESH_TEST_H

#include "multiregion_test.h"
int dg_mesh_cell_fetch_buffer_test(dg_mesh *mesh, int verbose);
int dg_mesh_parallel_test(dg_mesh *mesh, int verbose);
int dg_mesh_node_fetch_buffer_test(dg_mesh *mesh, int verbose);

#endif //DGOM_MR_MESH_TEST_H
