#include "SetTestMultiRegions.h"
#include "LibUtilities/GenUniformMesh.h"
#include "MultiRegBC2d_data.h"

#define N 1

MultiRegBC2d* setTriTestMesh(){
//    int N = 1;
    stdCell *shape = sc_create(N, TRIANGLE);
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    UnstructMesh *grid = ParallelUniformTriMesh_create(2, 2, -1, 1, -1, 1, 1, procid, nprocs);
    MultiReg2d *mesh = MultiReg2d_create(shape, grid);
    UnstructMesh_free(grid);

    extern int SFToV[NSURF][3];
    int **I_SFToV = IntMatrix_create(NSURF, 3);
    int k;
    for(k=0;k<NSURF;k++){
        I_SFToV[k][0] = SFToV[k][0];
        I_SFToV[k][1] = SFToV[k][1];
        I_SFToV[k][2] = SFToV[k][2];
    }
    MultiRegBC2d *surf = MultiRegBC2d_create(mesh, NSURF, I_SFToV);

    return surf;
}

MultiRegBC2d* setQuadTestMesh(){
//    int N = 1;
    stdCell *shape = sc_create(N, QUADRIL);
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    UnstructMesh *grid = ParallelUniformQuadMesh_create(2, 2, -1, 1, -1, 1, procid, nprocs);
    MultiReg2d *mesh = MultiReg2d_create(shape, grid);
    UnstructMesh_free(grid);

    extern int SFToV[NSURF][3];
    int **I_SFToV = IntMatrix_create(NSURF, 3);
    int k;
    for(k=0;k<NSURF;k++){
        I_SFToV[k][0] = SFToV[k][0];
        I_SFToV[k][1] = SFToV[k][1];
        I_SFToV[k][2] = SFToV[k][2];
    }
    MultiRegBC2d *surf = MultiRegBC2d_create(mesh, NSURF, I_SFToV);

    return surf;
}