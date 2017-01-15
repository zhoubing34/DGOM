//
// Created by li12242 on 16/12/16.
//

#include "pf_cellMean_test.h"
#include "PhysField/pf_cellMean.h"

int phys_cellMean_test(physField *phys, int verbose, char *message, char *filename){

    // local variable
    int fail = 0;

    multiReg *region = phys->region;
    parallMesh *mesh = phys->mesh;

    const int K = phys->grid->K;
    const int Np = phys->cell->Np;

    int k,i;
    // assignment
    int sk = 0;
    for(k=0;k<K;k++){
        for(i=0;i<Np;i++){
            phys->f_Q[sk++] = region->x[k][i];
            phys->f_Q[sk++] = region->y[k][i];
        }
    }

    double clockT1, clockT2;
    clockT1 = MPI_Wtime();
    pf_cellMean(phys);
    clockT2 = MPI_Wtime();

    if(verbose) {
        FILE *fp = CreateLog(filename, mesh->procid, mesh->nprocs);
        PrintVector2File(fp, "f_Q", phys->f_Q, K*Np*phys->Nfield);
        PrintVector2File(fp, "c_Q", phys->c_Q, K*phys->Nfield);
        fclose(fp);
    }

    return fail;
}
