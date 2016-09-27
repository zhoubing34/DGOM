#include "MultiRegionsTest.h"

int main(int argc, char **argv){

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    int N=3, Nvert=3;
    StdRegions2d *tri = StdTriEle_create(N);
    MultiReg2d *mesh;
    SetTestTriMesh(tri, mesh);

    printf("procid:%d, K = %d\n", mesh->procid, mesh->K);

    /* gen log filename */
    char casename[18] = "LoadBalanceTest";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

    /* write parEToV */
    PrintIntMatrix2File(fp, "EToV", mesh->EToV, mesh->K, Nvert);
    /* write GX */
    PrintMatrix2File(fp, "GX", mesh->GX, mesh->K, Nvert);
    /* write GY */
    PrintMatrix2File(fp, "GY", mesh->GY, mesh->K, Nvert);

    fclose(fp);
    MultiReg2d_free(mesh);
    StdRegions2d_free(tri);

    MPI_Finalize();
    return 0;
}

