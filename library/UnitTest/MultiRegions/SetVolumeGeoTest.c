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

    printf("init mesh\n");
    StdRegions2d *quad = GenStdQuadEle(N);
    MultiReg2d *mesh;
    SetTestQuadMesh(quad, mesh);

    printf("procid:%d, K = %d\n", mesh->procid, mesh->K);

    /* gen log filename */
    char casename[24] = "SetQuadVolumeGeoTest";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

    int n,m;
    /* write J */
    int sk=0;
    fprintf(fp, "J = \n");
    for(n=0;n<mesh->K;++n){
        for(m=0;m<quad->Np;++m){
            fprintf(fp, " %f, ", mesh->J[sk++]);
        }
        fprintf(fp, " \n");
    }

    /* write area */
    sk = 0;
    fprintf(fp, "Area = \n");
    for(n=0;n<mesh->K;++n){
        fprintf(fp, " %f, ", mesh->area[sk++]);
        fprintf(fp, " \n");
    }

    /* write J */
    sk=0;
    fprintf(fp, "vgeo = \n");
    for(n=0;n<mesh->K;++n){
        for(m=0;m<quad->Np;++m){
            fprintf(fp, "drdx = %f, drdy = %f, dsdx = %f, dsdy = %f\n",
                    mesh->vgeo[sk++],mesh->vgeo[sk++],
                    mesh->vgeo[sk++],mesh->vgeo[sk++]);
        }
        fprintf(fp, " \n");
    }
    fclose(fp);

    FreeMultiReg2d(mesh);
    FreeStdRegions2d(quad);
}

void TriTest(void){
    int N=3;

    printf("init mesh\n");
    StdRegions2d *tri = GenStdTriEle(N);
    MultiReg2d *mesh;
    SetTestTriMesh(tri, mesh);

    printf("procid:%d, K = %d\n", mesh->procid, mesh->K);

    /* gen log filename */
    char casename[24] = "SetTriVolumeGeoTest";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

    int n,m;
    /* write J */
    int sk=0;
    fprintf(fp, "J = \n");
    for(n=0;n<mesh->K;++n){
        for(m=0;m<tri->Np;++m){
            fprintf(fp, " %f, ", mesh->J[sk++]);
        }
        fprintf(fp, " \n");
    }

    /* write area */
    sk = 0;
    fprintf(fp, "Area = \n");
    for(n=0;n<mesh->K;++n){
        fprintf(fp, " %f, ", mesh->area[sk++]);
        fprintf(fp, " \n");
    }

    /* write J */
    sk=0;
    fprintf(fp, "vgeo = \n");
    for(n=0;n<mesh->K;++n){
        for(m=0;m<tri->Np;++m){
            fprintf(fp, "drdx = %f, drdy = %f, dsdx = %f, dsdy = %f\n",
                    mesh->vgeo[sk++],mesh->vgeo[sk++],
                    mesh->vgeo[sk++],mesh->vgeo[sk++]);
        }
        fprintf(fp, " \n");
    }

    fclose(fp);

    FreeMultiReg2d(mesh);
    FreeStdRegions2d(tri);
}