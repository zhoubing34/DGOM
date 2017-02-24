/**
 * @file
 * NetCDF变量与文件函数
 */

#include "nc_library.h"
#include "utility.h"

#define NCNAMELEN 1024

/**
 * @brief
 * create nc_dim variable
 * @param name name string
 * @param len length of dimension
 * @return
 * pointer to a new nc_dim
 */
nc_dim* nc_dim_create(const char *name, int len){
    nc_dim *dim = (nc_dim*) calloc(1, sizeof(nc_dim));

    /* allocate name */
    dim->name  = (char *) calloc(strlen(name), sizeof(char));
    /* assignment */
    strcpy(dim->name, name);
    if (len>0) {
        dim->len = len;
    }else if (len==0){
        dim->len = NC_UNLIMITED;
    }
    return dim;
}
/**
 * @brief free nc_dim variable
 * @param dim nc_dim variable
 */
static void nc_dim_free(nc_dim *dim){
    free(dim->name);
    free(dim);
}
/**
 * @brief print the infomation of nc_dim
 * @param dim nc_dim variable
 */
static void nc_dim_print(nc_dim *dim){
    printf("NetCDF Dimension\n");
    printf("  name:   %s\n", dim->name);
    printf("  length: %d\n", dim->len);
    printf("\n");
}
/**
 * @brief create nc_var
 * @param name variable name
 * @param ndim number of dimension
 * @param dim_vec_p pointer to the nc_dim vector
 * @param type variable type
 * @return var nc_var pointer
 */
nc_var* nc_var_create(const char *name, int ndim, nc_dim **dim_vec_p, int type){

    nc_var *var = (nc_var*) calloc(1, sizeof(nc_var));
    /* allocate name length */
    var->name  = (char *) calloc(strlen(name), sizeof(char));
    /* assignment */
    strcpy(var->name, name);
    var->ndim = ndim; /* number of dimensions */
    var->dim_vec_p = (nc_dim**) calloc(ndim, sizeof(nc_dim*) );
    int i;
    for(i=0;i<ndim;i++) {var->dim_vec_p[i] = dim_vec_p[i];}

    switch (type){
        case NC_INT: var->type = NC_INT; break;
        case NC_FLOAT: var->type = NC_FLOAT; break;
        case NC_DOUBLE: var->type = NC_DOUBLE; break;
        case NC_SHORT: var->type = NC_SHORT; break;
        case NC_STRING: var->type = NC_STRING; break;
        default: printf("%s(%d): Unknown type of NetCDF variable type: %d!\n",
                        __FILE__, __LINE__, type);
    }
    return var;
}
/**
 * @brief deallocate the nc_var
 * @param var nc_var variable
 */
static void nc_var_free(nc_var *var){
    free(var->name);
    free(var->dim_vec_p);
    free(var);
}
/**
 * @brief print the nc_var
 * @param var nc_var variable
 */
void nc_var_print(nc_var *var){
    printf("NetCDF Variable\n");
    printf("  name: %s\n", var->name);
    printf("  ndim: %d\n", var->ndim);
    printf("  type: %d\n", var->type);
    printf("  include dimensions\n");

    int i;
    for(i=0;i<var->ndim;i++){ nc_dim_print(var->dim_vec_p[i]); }
    printf("\n");
}
/**
 * @brief create new nc_file
 * @param name name of nc_file
 * @param procid process id
 * @param nprocs number of process
 * @param ndim number of nc_dim
 * @param dim_vec_p pointer to the nc_dim vector
 * @param nvar number of nc_var
 * @param var_vec_p pointer to the nc_var vector
 * @return file nc_file pointer
 */
nc_file* nc_file_create(const char *name, int procid, int nprocs,
                       int ndim, nc_dim **dim_vec_p,
                       int nvar, nc_var **var_vec_p){

    nc_file *file = (nc_file*) calloc(1, sizeof(nc_file));

    /* check name length */
    file->name  = (char *) calloc(strlen(name), sizeof(char));

    /* assignment */
    file->procid = procid; /* process id */
    file->nprocs = nprocs; /* process number */
    /* name */
    if(snprintf(file->name, NCNAMELEN, "%s%d-%d.nc", name, procid, nprocs) < 0 ){
        fprintf(stderr, "%s(%d): Error in creating NetCDF output file\n", __FILE__, __LINE__);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* dimension */
    file->ndim     = ndim;
    file->dim_vec_p = (nc_dim**) calloc(ndim, sizeof(nc_dim*) );
    int i;
    for(i=0;i<ndim;i++) {file->dim_vec_p[i] = dim_vec_p[i];}
    /* var */
    file->nvar     = nvar;
    file->var_vec_p = (nc_var**) calloc(nvar, sizeof(nc_var*) );
    for(i=0;i<nvar;i++) {file->var_vec_p[i] = var_vec_p[i];}

    return file;
}
/**
 * @brief print the info mation of nc_file
 * @param file nc_file pointer
 */
void nc_file_print(nc_file *file){
    printf("NetCDF Files\n");
    printf("  name:   %s\n", file->name);
    printf("  procid: %d\n", file->procid);
    printf("  nprocs: %d\n", file->nprocs);
    printf("  nvar:   %d\n", file->nvar);
    printf("  ndim:   %d\n", file->ndim);

    int i;
    for(i=0;i<file->nvar;i++){
        nc_var_print(file->var_vec_p[i]);
    }
    for(i=0;i<file->ndim;i++){
        nc_dim_print(file->dim_vec_p[i]);
    }
}
/**
 * @brief free the memory of nc_file and its included nc_dim and nc_var
 * @param file nc_file pointer
 */
void nc_file_free(nc_file *file){
    int i;
    free(file->name);
    for (i=0;i<file->ndim;i++)
        nc_dim_free(file->dim_vec_p[i]);

    for (i=0;i<file->nvar;i++)
        nc_var_free(file->var_vec_p[i]);

    free(file);
}

/**
 * @brief create the NetCDF file and define its dimensions and variables
 * @param file nc_file pointer
 */
void nc_file_init(nc_file *file){

    const int Ndim = file->ndim;
    const int Nvar = file->nvar;
    /* create output file */
    nc_error(ncmpi_create(MPI_COMM_SELF, file->name, NC_CLOBBER|NC_64BIT_OFFSET,
                          MPI_INFO_NULL, &file->id));
    /* define dimensions */
    int i,j;
    nc_dim *dim;
    for (i=0;i<Ndim;i++){
        dim = file->dim_vec_p[i]; /* pointer of dim */
        nc_error( ncmpi_def_dim(file->id, dim->name, dim->len, &dim->id) );
    }
    /* define variables */
    nc_var *var;
    for (i=0;i<Nvar; i++){
        var = file->var_vec_p[i]; /* pointer of var */
        int *dimids = (int *)calloc(var->ndim, sizeof(int));
        for(j=0;j<var->ndim;j++){ dimids[j] = var->dim_vec_p[j]->id; }

        nc_error(ncmpi_def_var(file->id, var->name, var->type,
                               var->ndim, dimids, &var->id));
        free(dimids);
    }
    /* finish definition */
    nc_error( ncmpi_enddef(file->id) );
}

/**
 * @brief close the NetCDF file
 * @param file nc_file pointer
 */
void nc_file_close(nc_file *file){
    nc_error( ncmpi_close(file->id) );

}