#include "MultiRegionsTest.h"

int TriCoordTest(void);
int QuadCoordTest(void);

int main(int argc, char **argv){
    int info;

    /* initialize MPI */
    MPI_Init(&argc, &argv);
    info = TriCoordTest();
    info = QuadCoordTest();

    MPI_Finalize();
    return info;
}

int QuadCoordTest(void){
    int info=0, N=1;

    printf("init tri mesh\n");
    StdRegions2d *quad = GenStdQuadEle(N);
    MultiReg2d *mesh;
    SetTestQuadMesh(quad, mesh);

    printf("procid:%d, K = %d\n", mesh->procid, mesh->K);

    /* gen log filename */
    char casename[24] = "SetQuadCoorNodeTest";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);
    PrintMatrix2File(fp, "gx", mesh->GX, mesh->K, mesh->stdcell->Nv);
    PrintMatrix2File(fp, "gy", mesh->GY, mesh->K, mesh->stdcell->Nv);

    fclose(fp);
    FreeMultiReg2d(mesh);
    FreeStdRegions2d(quad);
    return info;
}

int TriCoordTest(void){
    int info=0, N=1;

    printf("init tri mesh\n");
    StdRegions2d *tri = GenStdTriEle(N);
    MultiReg2d *mesh;
    SetTestTriMesh(tri, mesh);

    printf("procid:%d, K = %d\n", mesh->procid, mesh->K);

    /* gen log filename */
    char casename[24] = "SetTriCoorNodeTest";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);
    PrintMatrix2File(fp, "gx", mesh->GX, mesh->K, mesh->stdcell->Nv);
    PrintMatrix2File(fp, "gy", mesh->GY, mesh->K, mesh->stdcell->Nv);

    fclose(fp);
    FreeMultiReg2d(mesh);
    FreeStdRegions2d(tri);
    return info;
}