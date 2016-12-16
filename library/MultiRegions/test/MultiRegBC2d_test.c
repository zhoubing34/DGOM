//
// Created by li12242 on 16/12/15.
//

#include <MultiRegions/MultiRegBC/MultiRegBC2d.h>
#include "MultiRegBC2d_test.h"

int MultiTriRegions_MultiRegBC2d_test(MultiRegBC2d *surf, int verbose){
    int fail = 0;
    MultiReg2d *mesh = surf->mesh;
    int k;

    if(verbose) {
        char casename[34] = "MultiTriRegions_MultiRegBC2d_test";
        FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);
        PrintIntMatrix2File(fp, "EToBS", surf->EToBS, mesh->K, mesh->stdcell->Nv);
        fprintf(fp, "Nobc = %d\n", surf->Nobc);

        PrintIntVector2File(fp, "bcTypeList", surf->bcTypeList, surf->Nobc);
        for (k = 0; k < surf->Nobc; k++) {
            fprintf(fp, "surf[%d]->oblist->Nv = %d\n", k, surf->oblist[k]->Nv);
            PrintIntVector2File(fp, "BVToV", surf->oblist[k]->BVToV, surf->oblist[k]->Nv);
        }
        fclose(fp);
    }
    return fail;
}


int MultiQuadRegions_MultiRegBC2d_test(MultiRegBC2d *surf, int verbose){
    int fail = 0;
    MultiReg2d *mesh = surf->mesh;
    int k;

    if(verbose){
        char casename[35] = "MultiQuadRegions_MultiRegBC2d_test";
        FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);
        PrintIntMatrix2File(fp, "EToBS", surf->EToBS, mesh->K, mesh->stdcell->Nv);
        fprintf(fp, "Nobc = %d\n", surf->Nobc);
        PrintIntVector2File(fp, "bcTypeList", surf->bcTypeList, surf->Nobc);

        for(k=0;k<surf->Nobc;k++){
            fprintf(fp, "bc2d[%d]->oblist->Nv = %d\n", k, surf->oblist[k]->Nv);
            PrintIntVector2File(fp, "bc2d->oblist->BVToV", surf->oblist[k]->BVToV, surf->oblist[k]->Nv);
        }
        fclose(fp);
    }
    return fail;
}