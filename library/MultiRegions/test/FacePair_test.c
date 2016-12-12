//
// Created by li12242 on 12/11/16.
//

#include "FacePair_test.h"
#include "MulitRegions_test.h"
#include "SetTestMultiRegions.h"

int MultiTriRegions_FacePair_test(){
    int fail = 0;

    MultiReg2d* mesh = SetTriParallelMultiRegions();
    StdRegions2d *shape = mesh->stdcell;

    /* gen log filename */
    char casename[32] = "MultiTriRegions_FacePair_test";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

    /* write EToV */
    PrintIntMatrix2File(fp, "EToV", mesh->EToV, mesh->K, shape->Nv);
    /* write EToE */
    PrintIntMatrix2File(fp, "EToE", mesh->EToE, mesh->K, shape->Nv);
    /* write EToF */
    PrintIntMatrix2File(fp, "EToF", mesh->EToF, mesh->K, shape->Nv);
    /* write EToP */
    PrintIntMatrix2File(fp, "EToP", mesh->EToP, mesh->K, shape->Nv);

    fclose(fp);

    // free and finalize
    StdRegions2d_free(shape);
    MultiReg2d_free(mesh);
    return fail;
}

int MultiQuadRegions_FacePair_test(){
    int fail = 0;

    MultiReg2d* mesh = SetQuadParallelMultiRegions();
    StdRegions2d *shape = mesh->stdcell;

    /* gen log filename */
    char casename[32] = "MultiQuadRegions_FacePair_test";
    FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

    /* write EToV */
    PrintIntMatrix2File(fp, "EToV", mesh->EToV, mesh->K, shape->Nv);
    /* write EToE */
    PrintIntMatrix2File(fp, "EToE", mesh->EToE, mesh->K, shape->Nv);
    /* write EToF */
    PrintIntMatrix2File(fp, "EToF", mesh->EToF, mesh->K, shape->Nv);
    /* write EToP */
    PrintIntMatrix2File(fp, "EToP", mesh->EToP, mesh->K, shape->Nv);

    fclose(fp);

    // free and finalize
    StdRegions2d_free(shape);
    MultiReg2d_free(mesh);
    return fail;
}