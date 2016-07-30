#ifndef MULTIRESIONS_TEST_H
#define MULTIRESIONS_TEST_H

#include "MultiRegions/MultiRegions.h"
#include "LibUtilities/LibUtilities.h"

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
    int Test_i, Test_j, Kstart=0; \
    int procid, nprocs; \
    MPI_Comm_rank(MPI_COMM_WORLD, &procid); \
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); \
    if (nprocs<2) \
        TestMeshTri_Klocal[0] = 4; \
    for(Test_i=0;Test_i<procid;++Test_i){  \
        Kstart += TestMeshTri_Klocal[Test_i]; \
    } \
    int **parEToV = BuildIntMatrix(TestMeshTri_Klocal[procid], TestMeshTri_Nvert); \
    double *TestTri_VX = BuildVector(TestMeshTri_Nv); \
    double *TestTri_VY = BuildVector(TestMeshTri_Nv); \
    int Test_sk=0; \
    for(Test_i=0;Test_i<TestMeshTri_K;Test_i++){ \
        if(Test_i>=Kstart && Test_i<Kstart+TestMeshTri_Klocal[procid]) { \
            for (Test_j = 0; Test_j < TestMeshTri_Nvert; Test_j++) { \
                parEToV[Test_sk][Test_j] = TestMeshTri_EToV[Test_i][Test_j] - 1; } \
            Test_sk++; \
        } \
    } \
    mesh = GenMultiReg2d(tri, TestMeshTri_Klocal[procid], TestMeshTri_Nv, parEToV, TestMeshTri_VX, TestMeshTri_VY); \
}while(0) \

/* set quadrilateral test mesh */
#define SetTestQuadMesh(quad, mesh) \
do { \
    int Test_i, Test_j, Kstart=0; \
    int procid, nprocs; \
    MPI_Comm_rank(MPI_COMM_WORLD, &procid); \
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); \
    if (nprocs<2) \
        TestMeshQuad_Klocal[0] = 4; \
    for(Test_i=0;Test_i<procid;++Test_i){  \
        Kstart += TestMeshQuad_Klocal[Test_i]; \
    } \
    int **parEToV = BuildIntMatrix(TestMeshQuad_Klocal[procid], TestMeshQuad_Nvert); \
    int Test_sk=0; \
    for(Test_i=0;Test_i<TestMeshQuad_K;Test_i++){ \
        if(Test_i>=Kstart && Test_i<Kstart+TestMeshQuad_Klocal[procid]) { \
            for (Test_j = 0; Test_j < TestMeshQuad_Nvert; Test_j++) { \
                parEToV[Test_sk][Test_j] = TestMeshQuad_EToV[Test_i][Test_j] - 1; } \
            Test_sk++; \
        } \
    } \
    mesh = GenMultiReg2d(quad, TestMeshQuad_Klocal[procid], TestMeshQuad_Nv, parEToV, TestMeshQuad_VX, TestMeshQuad_VY); \
}while(0) \

#endif