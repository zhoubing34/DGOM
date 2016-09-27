#include "MultiRegionsTest.h"

void TriTest(void);
void QuadTest(void);
void PtintIntMat2File(FILE *fp, char *message, int **Mat, int row, int col);


int main(int argc, char **argv){
    /* initialize MPI */
    MPI_Init(&argc, &argv);

    TriTest();
    QuadTest();

    MPI_Finalize();
    return 0;
}

void QuadTest(void){
    int N=3, Nvert=4;

    StdRegions2d *quad = StdQuadEle_create(N);
    MultiReg2d *mesh;
    SetTestQuadMesh(quad, mesh);

    printf("procid:%d, K = %d\n", mesh->procid, mesh->K);

    /* gen log filename */
    char casename[24] = "SetQuadFacePairTest";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

    /* write EToV */
    PtintIntMat2File(fp, "EToV", mesh->EToV, mesh->K, Nvert);
    /* write EToE */
    PtintIntMat2File(fp, "EToE", mesh->EToE, mesh->K, Nvert);
    /* write EToF */
    PtintIntMat2File(fp, "EToF", mesh->EToF, mesh->K, Nvert);
    /* write EToP */
    PtintIntMat2File(fp, "EToP", mesh->EToP, mesh->K, Nvert);

    fclose(fp);

    MultiReg2d_free(mesh);
    StdRegions2d_free(quad);
}

void PtintIntMat2File(FILE *fp, char *message, int **Mat, int row, int col){
    fprintf(fp, "%s = \n", message);
    int n,m;
    for(n=0;n<row;++n){
        for(m=0;m<col;++m){
            fprintf(fp, " %d, ", Mat[n][m]);
        }
        fprintf(fp, " \n");
    }
}

void TriTest(void){
    int N=3, Nvert=3;

    printf("init mesh\n");
    StdRegions2d *tri = StdTriEle_create(N);
    MultiReg2d *mesh;
    SetTestTriMesh(tri, mesh);

    printf("procid:%d, K = %d\n", mesh->procid, mesh->K);

    /* gen log filename */
    char casename[24] = "SetTriFacePairTest";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

    /* write EToV */
    PtintIntMat2File(fp, "EToV", mesh->EToV, mesh->K, Nvert);
    /* write EToE */
    PtintIntMat2File(fp, "EToE", mesh->EToE, mesh->K, Nvert);
    /* write EToF */
    PtintIntMat2File(fp, "EToF", mesh->EToF, mesh->K, Nvert);
    /* write EToP */
    PtintIntMat2File(fp, "EToP", mesh->EToP, mesh->K, Nvert);

    fclose(fp);

    MultiReg2d_free(mesh);
    StdRegions2d_free(tri);
}