//
// Created by li12242 on 17/1/10.
//

#include "phys_limiter_test.h"
#include "PhysField/phys_limiter.h"

int phys_limiter_test(physField *phys, int verbose, char *message, char *filename){

    int fail = 0;

    multiReg *region = phys->region;
    const int K = phys->grid->K;
    const int Np = phys->cell->Np;

    int k,i;
    // assignment
    int sk = 0;
    for(k=0;k<K;k++){
        for(i=0;i<Np;i++){
            phys->f_Q[sk++] =      region->x[k][i];
            phys->f_Q[sk++] =     (region->y[k][i]);
        }
    }

    double clockT1, clockT2;
    clockT1 = MPI_Wtime();
    phys_slloc2d(phys, 2.0);
    clockT2 = MPI_Wtime();

    if(verbose) {
        FILE *fp = CreateLog(filename, phys->mesh->procid, phys->mesh->nprocs);
        PrintVector2File(fp, "f_Q", phys->f_Q, K*Np*phys->Nfield);
        PrintVector2File(fp, "c_Q", phys->c_Q, K*phys->Nfield);
        fclose(fp);
    }

    return fail;
}