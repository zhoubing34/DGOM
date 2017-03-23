//
// Created by li12242 on 17/2/18.
//

#include "dg_phys_init.h"

/***
 * @brief initialize the physical fields with initial condition file.
 * @param phys pointer to dg_phys structure;
 * @param casename name of case;
 * @note
 * The physical field is initialized with the vertex value
 */
void dg_phys_init_file2d(dg_phys *phys, char *casename){

    char init_filename[MAX_NAME_LENGTH];
    strcpy(init_filename, casename);
    strcat(init_filename, ".init");

    dg_grid *grid = dg_phys_grid(phys);

    /* open the init file */
    FILE *fp;
    dg_fopen(fp, init_filename, "Unable to open initial condition file");

    int Nvert, Nfield;
    fscanf(fp, "%d %d", &Nvert, &Nfield);
    if( Nvert != dg_grid_Nv(grid) ){
        fprintf(stderr, "%s (%d): Wrong number of vertex in initial file %s.\n",
                __FILE__, __LINE__, init_filename);
    }
    if( Nfield != dg_phys_Nfield(phys) ){
        fprintf(stderr, "%s (%d): Wrong physical field number in initial file %s.\n",
                __FILE__, __LINE__, init_filename);
    }

    /* read vertex initial values */
    register int k,n,fld;
    dg_real **init_value = matrix_real_create(Nvert, Nfield);
    for(n=0;n<Nvert;n++){
        int tmp; // node index
        fscanf(fp, "%d", &tmp);
        for(fld=0;fld<Nfield;fld++){
            fscanf(fp, "%lf", init_value[n]+fld);
        }
    }

    /* read element file */
    int **EToV = dg_grid_EToV(grid);

    /* assign to node fields */
    dg_cell *cell = dg_phys_cell(phys);
    const int Np = dg_cell_Np(cell);
    const int Nv = dg_cell_Nv(cell);
    const int K = dg_grid_K(grid);

    dg_real initial_vert_value[Nv*Nfield];
    dg_real initial_node_value[Np*Nfield];
    dg_real *f_Q = dg_phys_f_Q(phys);
    for(k=0;k<K;k++){
        for(n=0;n<Nv;n++){ // vertex initial values
            int vertID = EToV[k][n];
            for(fld=0;fld<Nfield;fld++){
                initial_vert_value[n*Nfield+fld] = init_value[vertID][fld];
            }
        }

        dg_cell_proj_vert2node(cell, Nfield, initial_vert_value, initial_node_value);
        int sk = (k*Np+n)*Nfield;
        for(n=0;n<Np*Nfield;n++){ // assign to node values
            f_Q[sk + n] = initial_node_value[n];
        }

    }

    matrix_real_free(init_value);
    matrix_int_free(EToV);
    return;
}