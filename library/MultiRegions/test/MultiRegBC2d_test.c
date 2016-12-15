//
// Created by li12242 on 16/12/15.
//

#include <MultiRegions/MultiRegBC/MultiRegBC2d.h>
#include <MultiRegions/MultiRegions.h>
#include "MultiRegBC2d_test.h"
#include "MultiRegBC2d_data.h"

int MultiTriRegions_MultiRegBC2d_test(MultiReg2d *mesh, int verbose){
    int fail = 0;

    extern int SFToV[NSURF][3];
    int **I_SFToV = IntMatrix_create(NSURF, 3);
    int k;
    for(k=0;k<NSURF;k++){
        I_SFToV[k][0] = SFToV[k][0];
        I_SFToV[k][1] = SFToV[k][1];
        I_SFToV[k][2] = SFToV[k][2];
    }
    MultiRegBC2d *bc2d = MultiRegBC2d_create(mesh, NSURF, I_SFToV);


    if(verbose) {
        char casename[34] = "MultiTriRegions_MultiRegBC2d_test";
        FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);
        PrintIntMatrix2File(fp, "EToBS", bc2d->EToBS, mesh->K, mesh->stdcell->Nv);
        fprintf(fp, "Nobc = %d\n", bc2d->Nobc);

        PrintIntVector2File(fp, "bcTypeList", bc2d->bcTypeList, bc2d->Nobc);
        for (k = 0; k < bc2d->Nobc; k++) {
            fprintf(fp, "bc2d[%d]->obc2d->Nv = %d\n", k, bc2d->obc2d[k]->Nv);
            PrintIntVector2File(fp, "BVToV", bc2d->obc2d[k]->BVToV, bc2d->obc2d[k]->Nv);
        }
        fclose(fp);
    }

    MultiRegBC2d_free(bc2d);
    IntMatrix_free(I_SFToV);
    return fail;
}


int MultiQuadRegions_MultiRegBC2d_test(MultiReg2d *mesh, int verbose){
    int fail = 0;

    extern int SFToV[NSURF][3];
    int **I_SFToV = IntMatrix_create(NSURF, 3);
    int k;
    for(k=0;k<NSURF;k++){
        I_SFToV[k][0] = SFToV[k][0];
        I_SFToV[k][1] = SFToV[k][1];
        I_SFToV[k][2] = SFToV[k][2];
    }
    MultiRegBC2d *bc2d = MultiRegBC2d_create(mesh, NSURF, I_SFToV);

    if(verbose){
        char casename[35] = "MultiQuadRegions_MultiRegBC2d_test";
        FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);
        PrintIntMatrix2File(fp, "EToBS", bc2d->EToBS, mesh->K, mesh->stdcell->Nv);
        fprintf(fp, "Nobc = %d\n", bc2d->Nobc);
        PrintIntVector2File(fp, "bcTypeList", bc2d->bcTypeList, bc2d->Nobc);

        for(k=0;k<bc2d->Nobc;k++){
            fprintf(fp, "bc2d[%d]->obc2d->Nv = %d\n", k, bc2d->obc2d[k]->Nv);
            PrintIntVector2File(fp, "bc2d->obc2d->BVToV", bc2d->obc2d[k]->BVToV, bc2d->obc2d[k]->Nv);
        }
        fclose(fp);
    }

    MultiRegBC2d_free(bc2d);
    IntMatrix_free(I_SFToV);
    return fail;
}