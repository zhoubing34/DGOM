#include "MulitRegions_test.h"
#include "SetTestMultiRegions.h"

int MultiTriRegions_NodePair_Test(){
    // local regions
    int fail = 0;
    MultiReg2d* mesh = SetTriParallelMultiRegions();
    StdRegions2d *shape = mesh->stdcell;


    // check
//    char casename[32] = "MultiTriRegions_NodePair_Test";
//    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);
//
//    /* write vmapM */
//    PrintIntVector2File(fp, "vmapM", mesh->vmapM, mesh->K*shape->Nfp * shape->Nfaces);
//    /* write vmapP */
//    PrintIntVector2File(fp, "vmapP", mesh->vmapP, mesh->K*shape->Nfp * shape->Nfaces);
//    fprintf(fp, "parNtotalout = %d\n", mesh->parNtotalout);
//    PrintIntVector2File(fp, "parmapOut", mesh->parmapOUT, mesh->parNtotalout);
//
//    fclose(fp);

    // free and finalize
    StdRegions2d_free(shape);
    MultiReg2d_free(mesh);

    return fail;
}

int MultiQuadRegions_NodePair_Test(){
    // local regions
    int fail = 0;
    MultiReg2d* mesh = SetQuadParallelMultiRegions();
    StdRegions2d *shape = mesh->stdcell;

    //check
    char casename[32] = "MultiQuadRegions_NodePair_Test";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

    /* write vmapM */
    PrintIntVector2File(fp, "vmapM", mesh->vmapM, mesh->K*shape->Nfp * shape->Nfaces);
    /* write vmapP */
    PrintIntVector2File(fp, "vmapP", mesh->vmapP, mesh->K*shape->Nfp * shape->Nfaces);
    fprintf(fp, "parNtotalout = %d\n", mesh->parNtotalout);
    PrintIntVector2File(fp, "parmapOut", mesh->parmapOUT, mesh->parNtotalout);

    fclose(fp);

    // free and finalize
    StdRegions2d_free(shape);
    MultiReg2d_free(mesh);
    return fail;
}