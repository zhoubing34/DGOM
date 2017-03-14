//
// Created by li12242 on 17/2/18.
//

#include "dg_phys_init.h"

static int **read_EToV_file(dg_phys *phys, char *casename);

/***
 * @brief initialize the physical fields with initial condition file.
 * @param phys physical structure
 * @param casename name of computational case
 * @note
 * The EToV in phys->grid structure is different from the original EToV
 * variable. Therefore, a new EToV is read from the input element file
 * and assigned in the function.
 */
void dg_phys_init_file2d(dg_phys *phys, char *casename){

    char init_filename[MAX_NAME_LENGTH];
    strcpy(init_filename, casename);
    strcat(init_filename, ".init");

    /* open the init file */
    FILE *fp;
    dg_fopen(fp, init_filename, "Unable to open initial condition file");

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
    int **EToV = read_EToV_file(phys, casename);

    /* assign to node fields */
    dg_cell *cell = phys->cell;
    const int Np = dg_cell_Np(cell);
    const int Nv = dg_cell_Nv(cell);
    const int K = dg_grid_K(phys->grid);

    dg_real initial_vertex_value[Nv], initial_node_value[Np];
    dg_real *f_Q = phys->f_Q;
    for(k=0;k<K;k++){
        for(fld=0;fld<Nfield;fld++){
            for(n=0;n<Nv;n++){ // vertex initial values
                initial_vertex_value[n] = init_value[ EToV[k][n] ][fld];
            }
            /* map from vertex to nodes */
            dg_cell_proj_vert2node(cell, initial_vertex_value, initial_node_value);
            for(n=0;n<Np;n++){ // assign to node values
                f_Q[(k*Np+n)*Nfield + fld] = initial_node_value[n];
            }
        }
    }

    matrix_real_free(init_value);
    matrix_int_free(EToV);
    return;
}

static int **read_EToV_file(dg_phys *phys, char *casename){

    char element_file[MAX_NAME_LENGTH];
    strcpy(element_file, casename);
    strcat(element_file, ".ele");
    FILE *fp;
    if( (fp = fopen(element_file, "r")) == NULL ){
        fprintf(stderr, "%s (%d)\n"
                        "Unable to open element file %s.\n",
                __FUNCTION__,__LINE__,element_file);
    }
    int Nv, K, temp;
    // read cell number and cell vertex number
    fscanf(fp, "%d %d %d\n", &K, &Nv, &temp);
    // check element vertex
    if(dg_cell_Nv(phys->cell) !=  Nv){
        fprintf(stderr, "%s (%d)\n"
                "The input element type is not correct!\n", __FILE__, __LINE__);
        exit(-1);
    }
    int **EToV = matrix_int_create(K, Nv);
    int n,k;
    for(k=0;k<K;k++){
        fscanf(fp, "%d", &temp); //read index
        for(n=0;n<Nv;n++){
            fscanf(fp, "%d", EToV[0]+k*Nv+n);
            EToV[k][n] -= 1; // change index start from 0 (C style)
        }
        fscanf(fp, "%d", &temp); //read region id
    }
    fclose(fp);
    return EToV;
}