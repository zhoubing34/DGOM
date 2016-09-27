#ifndef MULTIRESIONS_TEST_H
#define MULTIRESIONS_TEST_H

#include "MultiRegions/MultiRegions.h"
#include "LibUtilities/LibUtilities.h"
#include "LibUtilities/GenUniformMesh.h"

int    TestMeshTri_Nv         = 6;
int    TestMeshTri_Nvert      = 3;
int    TestMeshTri_EToV[4][3] = {{1,4,2},{2,4,5},{2,5,6},{2,6,3}};
double TestMeshTri_VX[6]      = {-1, 0, 1, -1, 0, 1};
double TestMeshTri_VY[6]      = { 1, 1, 1,  0, 0, 0};
int    TestMeshTri_K          = 4;
int    TestMeshTri_Klocal[2]  = {2, 2};

int    TestMeshQuad_Nv        = 9;
int    TestMeshQuad_Nvert     = 4;
int    TestMeshQuad_EToV[4][4]= {{1,4,5,2},{2,5,6,3},{4,7,8,5},{5,8,9,6}};
double TestMeshQuad_VX[9]     = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
double TestMeshQuad_VY[9]     = {1, 1, 1, 0, 0, 0, -1, -1, -1};
int    TestMeshQuad_K          = 4;
int    TestMeshQuad_Klocal[2]  = {2, 2};

/* set triangle test mesh */
#define SetTestTriMesh(tri, mesh) \
do { \
    int procid, nprocs; \
    MPI_Comm_rank(MPI_COMM_WORLD, &procid); \
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); \
    UnstructMesh *grid; \
    if (nprocs<2){ \
        grid = GenUniformTriMesh(2, 2, -1, 1, -1, 1, 1); \
    }else{ \
        grid = GenParallelUniformTriMesh(2, 2, -1, 1, -1, 1, 1, procid, nprocs); \
    } \
    mesh = GenMultiReg2d(tri, grid); \
    DestroyUnstructMesh(grid); \
}while(0) \

/* set quadrilateral test mesh */
#define SetTestQuadMesh(quad, mesh) \
do { \
    int procid, nprocs; \
    MPI_Comm_rank(MPI_COMM_WORLD, &procid); \
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); \
    UnstructMesh *grid; \
    if (nprocs<2){ \
        grid = GenUniformQuadMesh(2, 2, -1, 1, -1, 1); \
    }else{ \
        grid = GenParallelUniformQuadMesh(2, 2, -1, 1, -1, 1, procid, nprocs); \
    } \
    mesh = GenMultiReg2d(quad, grid); \
    DestroyUnstructMesh(grid); \
}while(0) \

#endif