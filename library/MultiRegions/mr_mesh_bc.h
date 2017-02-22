//
// Created by li12242 on 12/21/16.
//

#ifndef DGOM_MR_MESH_ADDBOUNDARY_H
#define DGOM_MR_MESH_ADDBOUNDARY_H

#include "mr_mesh.h"
void mr_mesh_addBoundary2d(parallMesh *mesh, int Nsurf, int **SFToV);
void mr_mesh_deleteBoundary2d(parallMesh *mesh);
void mr_mesh_read_bcfile2d(parallMesh *mesh, char *casename);

#endif //DGOM_MR_MESH_ADDBOUNDARY_H