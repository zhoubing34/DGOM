//
// Created by li12242 on 17/2/18.
//

#include "pf_init_file.h"

/***
 * @brief initialize the physical fields with initial condition file.
 * @param phys physical structure
 * @param casename name of computational case
 */
void pf_init_file2d(physField *phys, char *casename){
    char init_filename[MAX_NAME_LENGTH];
    strcpy(init_filename, casename);
    strcat(init_filename, ".init");
    FILE *fp;
    /* open the init file */
    if( (fp = fopen(init_filename, "r")) == NULL ){
        fprintf(stderr, "pf_init_file2d (%s): %d\n"
                        "Unable to open initial condition file %s.\n",
                __FILE__,__LINE__,init_filename);
        exit(-1);
    }

    int Nvert, Nfield;
    fscanf(fp, "%d %d", &Nvert, &Nfield);
    if( Nvert != phys->grid->Nv ){
        fprintf(stderr, "%s (%d): Wrong number of vertex in initial file %s.\n",
                __FILE__, __LINE__, init_filename);
    }
    if( Nfield != phys->Nfield ){
        fprintf(stderr, "%s (%d): Wrong physical field number in initial file %s.\n",
                __FILE__, __LINE__, init_filename);
    }

    /* read vertex initial values */
    register int n,fld,k;
    int tmp;
    dg_real **intval = matrix_real_create(Nfield, Nvert);
    for(n=0;n<Nvert;n++){
        fscanf(fp, "%d", &tmp);
        for(fld=0;fld<Nfield;fld++){
            fscanf(fp, "%lf", intval[fld]+n);
        }
    }

    /* assign to node fields */
    stdCell *cell = phys->cell;
    const int Nv = phys->cell->Nv;
    const int K = phys->grid->K;
    const int Np = phys->cell->Np;
    int **EToV = phys->grid->EToV;
    dg_real floc_v[Nv], floc[Np];
    dg_real *f_Q = phys->f_Q;
    for(k=0;k<K;k++){
        for(fld=0;fld<Nfield;fld++){
            for(n=0;n<Nv;n++){ // vertex initial values
                floc_v[n] = intval[fld][ EToV[k][n] ];
            }
            /* map from vertex to nodes */
            sc_vertProj(cell, floc_v, floc);
            for(n=0;n<Np;n++){ // assign to node values
                f_Q[(k*Np+n)*Nfield + fld] = floc[n];
            }
        }
    }
    matrix_real_free(intval);
    return;
}