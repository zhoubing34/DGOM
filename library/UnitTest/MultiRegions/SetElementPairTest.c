//
// Created by li12242 on 16/7/30.
//

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

    printf("init quad mesh\n");
    StdRegions2d *quad = StdQuadEle_create(N);
    MultiReg2d *mesh;
    SetTestQuadMesh(quad, mesh);

    printf("procid:%d, K = %d\n", mesh->procid, mesh->K);

    /* gen log filename */
    char casename[24] = "SetQuadElementPairTest";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

    fprintf(fp, "parEtotalout = %d\n", mesh->parEtotalout);
    PrintIntVector2File(fp, "elemapOUT", mesh->elemapOut, mesh->parEtotalout);

    fclose(fp);
    MultiReg2d_free(mesh);
    StdRegions2d_free(quad);
}

void TriTest(void){
    int N=3;

    printf("init tri mesh\n");
    StdRegions2d *tri = StdTriEle_create(N);
    MultiReg2d *mesh;
    SetTestTriMesh(tri, mesh);

    printf("procid:%d, K = %d\n", mesh->procid, mesh->K);

    /* gen log filename */
    char casename[24] = "SetTriElementPairTest";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

    fprintf(fp, "parEtotalout = %d\n", mesh->parEtotalout);
    PrintIntVector2File(fp, "elemapOUT", mesh->elemapOut, mesh->parEtotalout);

    fclose(fp);

    MultiReg2d_free(mesh);
    StdRegions2d_free(tri);
}

