#ifndef UNSTRCUTMESH_H
#define UNSTRCUTMESH_H

#include "StdRegions/StdRegions.h"

typedef struct{
    int dim;                    // Dimensions
    int nv;                     // Num of vertex
    int ne;                     // Num of element
    ElementType eletype;        // type of element
    int **EToV;                 // vertex list in elements, the index start from 0
    double *vx, *vy, *vz;       // coordinate of vertex
} UnstructMesh;

void UnstructMesh_free(UnstructMesh *grid);
UnstructMesh* UnstructMesh_create(int Dim, int Ne, int Nv, ElementType eletype);

#endif // UNSTRCUTMESH_H
