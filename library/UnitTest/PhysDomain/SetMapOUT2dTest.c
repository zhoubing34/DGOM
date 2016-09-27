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
    StdRegions2d *quad = StdQuadEle_create(N);
    MultiReg2d *mesh;
    SetTestQuadMesh(quad, mesh);
    PhysDomain2d *phys = GenPhysDomain2d(mesh, Nfields);

    printf("procid:%d, K = %d\n", mesh->procid, mesh->K);
    /* gen log filename */
    char casename[24] = "QuadSetMapOut2dTest";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

    /* write mesh->parmapOUT */
    PrintIntVector2File(fp, "quad->parmapOUT", mesh->parmapOUT, mesh->parNtotalout);
    /* write phys->parmapOUT */
    PrintIntVector2File(fp, "QuadPhys->parmapOUT", phys->parmapOUT, phys->parNtotalout);

    fclose(fp);

    MultiReg2d_free(mesh);
    StdRegions2d_free(quad);
}

void TriDomainTest(void){
    int N=3, Nfields=3;

    printf("init mesh\n");
    StdRegions2d *tri = StdTriEle_create(N);
    MultiReg2d *mesh;
    SetTestTriMesh(tri, mesh);
    PhysDomain2d *phys = GenPhysDomain2d(mesh, Nfields);

    printf("procid:%d, K = %d\n", mesh->procid, mesh->K);

    /* gen log filename */
    char casename[24] = "TriSetMapOut2dTest";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

    /* write mesh->parmapOUT */
    PrintIntVector2File(fp, "quad->parmapOUT", mesh->parmapOUT, mesh->parNtotalout);
    /* write phys->parmapOUT */
    PrintIntVector2File(fp, "QuadPhys->parmapOUT", phys->parmapOUT, phys->parNtotalout);

    fclose(fp);

    MultiReg2d_free(mesh);
    StdRegions2d_free(tri);
}



