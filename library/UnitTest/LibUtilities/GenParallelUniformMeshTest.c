#include "LibUtiltiesTest.h"

int main(int argc, char **argv){

    int procid, nprocs;

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* gen log filename */
    char casename[28] = "GenParallelUniformMeshTest";
    FILE *fp = CreateLog(casename, procid, nprocs);


    double xmin = -1, xmax = 1;
    double ymin = -1, ymax = 1;

//    int **parEToV;
//    double *VX, *VY;
//    int Nv, Ne;
    UnstructMesh *grid;

    grid = ParallelUniformTriMesh_create(2, 2, xmin, xmax, ymin, ymax, 1, procid, nprocs);

    PrintIntMatrix2File(fp, "EToV for tri", grid->EToV, grid->ne, 3);
    PrintVector2File(fp, "VX for tri", grid->vx, grid->nv);
    PrintVector2File(fp, "VY for tri", grid->vy, grid->nv);

    UnstructMesh_free(grid);

    grid = ParallelUniformQuadMesh_create(2, 2, xmin, xmax, ymin, ymax, procid, nprocs);

    PrintIntMatrix2File(fp, "EToV for quad", grid->EToV, grid->ne, 4);
    PrintVector2File(fp, "VX for quad", grid->vx, grid->nv);
    PrintVector2File(fp, "VY for quad", grid->vy, grid->nv);

    UnstructMesh_free(grid);

    MPI_Finalize();

    return 0;
}

