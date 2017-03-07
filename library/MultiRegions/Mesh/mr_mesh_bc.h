//
// Created by li12242 on 12/21/16.
//

#ifndef DGOM_MR_MESH_ADDBOUNDARY_H
#define DGOM_MR_MESH_ADDBOUNDARY_H

#include "mr_mesh.h"
void mr_mesh_add_bc2d(dg_mesh *mesh, int Nsurf, int **SFToV);
void mr_mesh_del_bc2d(dg_mesh *mesh);
void mr_mesh_read_bcfile2d(dg_mesh *mesh, char *casename);

#endif //DGOM_MR_MESH_ADDBOUNDARY_H
