#include "MultiRegionsTest.h"

int main(int argc, char **argv){

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    int i,j, N=3, Nvert=3;
    int K=4, Nv = 6;
    int Np = (N+1)*(N+2)/2;
    int info, procid, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int EToVt[4][3] = {{1,4,2},{2,6,3},
                       {2,4,5},{2,5,6}};
    double VXt[6] = {-1, 0, 1, -1, 0, 1};
    double VYt[6] = {1, 1, 1, 0, 0, 0};

    int Kprocs[2] = {2,2};
    int Klocal[2] = {2,2};
    int Kstart=0;

    /* start index of element */
    for(i=0;i<procid;++i)
        Kstart += Kprocs[i];

    printf("procid:%d, Kstart = %d\n", procid, Kstart);
    printf("procid:%d, Kstart+Klocal = %d\n", procid, Kstart+Klocal[procid]);

    int **EToV = BuildIntMatrix(Kprocs[procid], Nvert);
    double *VX = BuildVector(Nv);
    double *VY = BuildVector(Nv);

    for(i=0;i<Nv;i++){
        VX[i] = VXt[i];
        VY[i] = VYt[i];
    }

    int sk=0;
    for(i=0;i<K;i++){
        if(i>=Kstart && i<Kstart+Klocal[procid]) {
            for (j = 0; j < Nvert; j++) {
                EToV[sk][j] = EToVt[i][j] - 1;
            }
            sk++;
        }
    }

    printf("init mesh\n");
    StdRegions2d *tri = GenStdTriEle(N);
    MultiReg2d *mesh = GenMultiReg2d(tri, Klocal[procid], Nv, EToV, VX, VY);

    printf("procid:%d, K = %d\n", mesh->procid, mesh->K);

    /* gen log filename */
    char casename[18] = "SetFacePairTest";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

    int n,m;
    /* write EToV */
    fprintf(fp, "EToV = \n");
    for(n=0;n<mesh->K;++n){
        for(m=0;m<Nvert;++m){
            fprintf(fp, " %d, ", mesh->EToV[n][m]);
        }
        fprintf(fp, " \n");
    }
    /* write EToE */
    fprintf(fp, "EToE = \n");
    for(n=0;n<mesh->K;++n){
        for(m=0;m<Nvert;++m){
            fprintf(fp, " %d, ", mesh->EToE[n][m]);
        }
        fprintf(fp, " \n");
    }
    /* write EToF */
    fprintf(fp, "EToF = \n");
    for(n=0;n<mesh->K;++n){
        for(m=0;m<Nvert;++m){
            fprintf(fp, " %d, ", mesh->EToF[n][m]);
        }
        fprintf(fp, " \n");
    }

    /* write EToP */
    fprintf(fp, "EToP = \n");
    for(n=0;n<mesh->K;++n){
        for(m=0;m<Nvert;++m){
            fprintf(fp, " %d, ", mesh->EToP[n][m]);
        }
        fprintf(fp, " \n");
    }

    fclose(fp);

    FreeMultiReg2d(mesh);
    FreeStdRegions2d(tri);

    MPI_Finalize();
    return 0;
}