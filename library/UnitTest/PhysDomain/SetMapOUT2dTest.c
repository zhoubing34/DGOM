#include "PhysDomainTest.h"
#include "../MultiRegions/MultiRegionsTest.h"

int main(int argc, char **argv){

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    void TriDomainTest(void);
    TriDomainTest();

    void QuadDomainTest(void);
    QuadDomainTest();

    MPI_Finalize();
    return 0;
}


void QuadDomainTest(void){
    int N=3, Nfields=3;

    printf("init mesh\n");
    StdRegions2d *quad = GenStdQuadEle(N);
    MultiReg2d *mesh;
    SetTestQuadMesh(quad, mesh);
    PhysDomain2d *phys = GetPhysDomain2d(mesh, Nfields);

    printf("procid:%d, K = %d\n", mesh->procid, mesh->K);
    /* gen log filename */
    char casename[24] = "QuadSetMapOut2dTest";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

    /* write mesh->parmapOUT */
    PrintIntVector2File(fp, "quad->parmapOUT", mesh->parmapOUT, mesh->parNtotalout);
    /* write phys->parmapOUT */
    PrintIntVector2File(fp, "QuadPhys->parmapOUT", phys->parmapOUT, phys->parNtotalout);

    fclose(fp);

    FreeMultiReg2d(mesh);
    FreeStdRegions2d(quad);
}

void TriDomainTest(void){
    int N=3, Nfields=3;

    printf("init mesh\n");
    StdRegions2d *tri = GenStdTriEle(N);
    MultiReg2d *mesh;
    SetTestTriMesh(tri, mesh);
    PhysDomain2d *phys = GetPhysDomain2d(mesh, Nfields);

    printf("procid:%d, K = %d\n", mesh->procid, mesh->K);

    /* gen log filename */
    char casename[24] = "TriSetMapOut2dTest";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

    /* write mesh->parmapOUT */
    PrintIntVector2File(fp, "quad->parmapOUT", mesh->parmapOUT, mesh->parNtotalout);
    /* write phys->parmapOUT */
    PrintIntVector2File(fp, "QuadPhys->parmapOUT", phys->parmapOUT, phys->parNtotalout);

    fclose(fp);

    FreeMultiReg2d(mesh);
    FreeStdRegions2d(tri);
}



