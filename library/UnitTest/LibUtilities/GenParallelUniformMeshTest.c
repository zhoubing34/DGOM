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

    int **parEToV;
    double *VX, *VY;
    int Nv, Ne;

    GenParallelUniformTriMesh(2, 2,
                              xmin, xmax, ymin, ymax, 1,
                              procid, nprocs,
                              &Ne, &Nv, &parEToV, &VX, &VY);

    PrintIntMatrix2File(fp, "EToV for tri", parEToV, Ne, 3);
    PrintVector2File(fp, "VX for tri", VX, Nv);
    PrintVector2File(fp, "VY for tri", VY, Nv);

    DestroyIntMatrix(parEToV);
    DestroyVector(VX);
    DestroyVector(VY);

    GenParallelUniformQuadMesh(2, 2, xmin, xmax, ymin, ymax, procid, nprocs, &Ne, &Nv, &parEToV, &VX, &VY);

    PrintIntMatrix2File(fp, "EToV for quad", parEToV, Ne, 4);
    PrintVector2File(fp, "VX for quad", VX, Nv);
    PrintVector2File(fp, "VY for quad", VY, Nv);

    DestroyIntMatrix(parEToV);
    DestroyVector(VX);
    DestroyVector(VY);

    return 0;
}

