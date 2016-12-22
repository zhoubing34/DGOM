//
// Created by li12242 on 16/12/16.
//

#include "CellFetch_test.h"

int cellFetch_test(PhysDomain2d *phys, int verbose, char *message, char *filename){

    // local variable
    int fail = 0;

    MultiReg2d *mesh = phys->mesh;
    MultiRegBC2d *surf = phys->surf;
    stdCell *shape = mesh->stdcell;

    int K = mesh->K;
    int Np = shape->Np;
    int Nfp = shape->Nfp;
    int Nfaces = shape->Nfaces;
    int nprocs = mesh->nprocs;

    int Nnode = Nfaces*Nfp*K;
    double *xM = Vector_create(K);
    double *xP = Vector_create(K);
    double *yM = Vector_create(K);
    double *yP = Vector_create(K);



}
