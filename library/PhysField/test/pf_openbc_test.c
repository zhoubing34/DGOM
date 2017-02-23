//
// Created by li12242 on 17/2/23.
//

#include <PhysField/pf_phys.h>
#include "pf_openbc_test.h"

int pf_openbc_test(physField *phys, int verbose){
    const int procid = phys->mesh->procid;
    const int nprocs = phys->mesh->nprocs;
    const int Nfield = phys->Nfield;
    const int K = phys->grid->K;
    const int Np = phys->cell->Np;
    int fail=0;

    int k,n,fld,sk=0;
    for(k=0;k<K;k++){
        for(n=0;n<Np;n++){
            for(fld=0;fld<Nfield;fld++){
                phys->f_ext[sk++] = 0.0;
            }
        }
    }

    pf_set_openbc(phys, 0, time_interp_linear);
    pf_set_openbc(phys, 0.3, time_interp_linear);
    pf_set_openbc(phys, 1, time_interp_linear);

    if(verbose){
        FILE *fp = create_log(__FUNCTION__, procid, nprocs);
        PrintVector2File(fp, "f_ext=", phys->f_ext, Nfield*Np*K);
        fclose(fp);
    }

    if(!procid){
        if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    }
    return fail;
}