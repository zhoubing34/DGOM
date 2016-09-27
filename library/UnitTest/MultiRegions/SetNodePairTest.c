#include "MultiRegionsTest.h"

void TriTest(void);
void QuadTest(void);

int main(int argc, char **argv){

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    TriTest();
    QuadTest();

    MPI_Finalize();
    return 0;
}

void QuadTest(void){
    int N=3;

    printf("init tri mesh\n");
    StdRegions2d *quad = StdQuadEle_create(N);
    MultiReg2d *mesh;
    SetTestQuadMesh(quad, mesh);

    printf("procid:%d, K = %d\n", mesh->procid, mesh->K);

    /* gen log filename */
    char casename[24] = "SetQuadNodePairTest";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

    /* write vmapM */
    PrintIntVector2File(fp, "vmapM", mesh->vmapM, mesh->K*quad->Nfp * quad->Nfaces);
    /* write vmapP */
    PrintIntVector2File(fp, "vmapP", mesh->vmapP, mesh->K*quad->Nfp * quad->Nfaces);
    fprintf(fp, "parNtotalout = %d\n", mesh->parNtotalout);
    PrintIntVector2File(fp, "parmapOut", mesh->parmapOUT, mesh->parNtotalout);

    fclose(fp);
    MultiReg2d_free(mesh);
    StdRegions2d_free(quad);
}

void TriTest(void){
    int N=3;

    printf("init quad mesh\n");
    StdRegions2d *tri = StdTriEle_create(N);
    MultiReg2d *mesh;
    SetTestTriMesh(tri, mesh);

    printf("procid:%d, K = %d\n", mesh->procid, mesh->K);

    /* gen log filename */
    char casename[24] = "SetTriNodePairTest";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

    /* write vmapM */
    PrintIntVector2File(fp, "vmapM", mesh->vmapM, mesh->K*tri->Nfp * tri->Nfaces);
    /* write vmapP */
    PrintIntVector2File(fp, "vmapP", mesh->vmapP, mesh->K*tri->Nfp * tri->Nfaces);
    fprintf(fp, "parNtotalout = %d\n", mesh->parNtotalout);
    PrintIntVector2File(fp, "parmapOut", mesh->parmapOUT, mesh->parNtotalout);
    fclose(fp);

    MultiReg2d_free(mesh);
    StdRegions2d_free(tri);
}