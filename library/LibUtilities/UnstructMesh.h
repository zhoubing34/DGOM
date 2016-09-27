#ifndef UNSTRCUTMESH_H
#define UNSTRCUTMESH_H

#include "LibUtilities/LibUtilities.h"

typedef enum {
    TRIANGLE, // triangle
    QUADRIL,  // quadrilateral
    TETRA,    // tetrahedron
    TRIPRISM, // triangle prism
    HEXA      // hexahedron
} EleType;

typedef struct{
    int dim;                    // Dimensions
    int nv;                     // Num of vertex
    int ne;                     // Num of element
    EleType eletype;            // type of each element
    int **EToV;                 // vertex list in elements, the index start from 0
    double *vx, *vy, *vz;       // coordinate of vertex
    char *name;   // name string
} UnstructMesh;


void DestroyUnstructMesh(UnstructMesh *grid);
UnstructMesh* CreateUnstructMesh(
        int Dim, int Ne, int Nv, EleType eletype);

#endif // UNSTRCUTMESH_H
