#ifndef UNSTRCUTMESH_H
#define UNSTRCUTMESH_H

#include "LibUtilities/LibUtilities.h"

typedef enum {
    TRIANGLE, // triangle
    QUADRIL,  // quadrilateral
    TETRA,    // tetrahedron
    TRIPRISM, // triangle prism
    HEXA      // hexahedron
} ElementType;

typedef struct{
    int dim;                    // Dimensions
    int nv;                     // Num of vertex
    int ne;                     // Num of element
    ElementType eletype;            // type of each element
    int **EToV;                 // vertex list in elements, the index start from 0
    double *vx, *vy, *vz;       // coordinate of vertex
    char *name;   // name string
} UnstructMesh;


void UnstructMesh_free(UnstructMesh *grid);
UnstructMesh* UnstructMesh_create(int Dim, int Ne, int Nv, ElementType eletype);

#endif // UNSTRCUTMESH_H
