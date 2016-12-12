#include "SetTestMultiRegions.h"
#include "LibUtilities/GenUniformMesh.h"

#define N 1

MultiReg2d* SetTriParallelMultiRegions(){
//    int N = 1;
    StdRegions2d *shape = StdRegions2d_create(N, TRIANGLE);
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    UnstructMesh *grid = ParallelUniformTriMesh_create(2, 2, -1, 1, -1, 1, 1, procid, nprocs);
    MultiReg2d *mesh = MultiReg2d_create(shape, grid);
    UnstructMesh_free(grid);
    return mesh;
}

MultiReg2d* SetQuadParallelMultiRegions(){
//    int N = 1;
    StdRegions2d *shape = StdRegions2d_create(N, QUADRIL);
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    UnstructMesh *grid = ParallelUniformQuadMesh_create(2, 2, -1, 1, -1, 1, procid, nprocs);
    MultiReg2d *mesh = MultiReg2d_create(shape, grid);
    UnstructMesh_free(grid);
    return mesh;
}