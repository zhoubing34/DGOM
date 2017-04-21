/**
 * @file
 * NetCDF变量与文件函数
 */

#include "nc_library.h"
#include "utility.h"
#include "unit_test.h"

#define NC_MAX_NAME_LEN 1024

/**
 * @brief
 * create nc_dim variable
 * @param name name string;
 * @param len length of dimension;
 * @return
 * pointer to a new nc_dim
 */
NC_Dim* nc_dim_create(const char *name, int len){
    NC_Dim *dim = (NC_Dim*) calloc(1, sizeof(NC_Dim));

    /* allocate name */
    dim->name  = (char *) calloc(strlen(name), sizeof(char));
    /* assignment */
    strcpy(dim->name, name);
    if (len>0) { dim->len = len; }
    else if (len==0){ dim->len = NC_UNLIMITED; }
    return dim;
}
/**
 * @brief free nc_dim variable
 * @param dim nc_dim variable
 */
static void nc_dim_free(NC_Dim *dim){
    free(dim->name);
    free(dim);
}
/**
 * @brief print the infomation of nc_dim
 * @param dim nc_dim variable
 */
static void nc_dim_print(NC_Dim *dim){
    printf(HEADEND "NetCDF Dimension\n");
    printf(HEADLINE "  name:   %s\n", dim->name);
    printf(HEADLINE "  length: %d\n", dim->len);
}
/**
 * @brief create nc_var
 * @param name variable name
 * @param ndim number of dimension
 * @param dim_vec_p pointer to the nc_dim vector
 * @param type variable type
 * @return var nc_var pointer
 */
NC_Var* nc_var_create(const char *name, int ndim, NC_Dim **dim_vec_p, int type){

    NC_Var *var = (NC_Var*) calloc(1, sizeof(NC_Var));
    /* allocate name length */
    var->name  = (char *) calloc(strlen(name), sizeof(char));
    /* assignment */
    strcpy(var->name, name);
    var->ndim = ndim; /* number of dimensions */
    var->dim_vec_p = (NC_Dim**) calloc(ndim, sizeof(NC_Dim*) );
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
static void nc_var_free(NC_Var *var){
    free(var->name);
    free(var->dim_vec_p);
    free(var);
}
/**
 * @brief print the nc_var
 * @param var nc_var variable
 */
void nc_var_print(NC_Var *var){
    printf(HEADEND "NetCDF Variable\n");
    printf(HEADLINE "  name: %s\n", var->name);
    printf(HEADLINE "  Ndim: %d\n", var->ndim);
    printf(HEADLINE "  type: %d\n", var->type);
    printf(HEADEND "Including dimensions:\n");
    int i;
    for(i=0;i<var->ndim;i++){
        printf(HEADLINE "  name: %s\n", var->dim_vec_p[i]->name);
    }
}
/**
 * @brief create pointer to a new nc_file structure.
 * @param name name of NC_File;
 * @param ndim number of NC_Dim;
 * @param dim_vec_p pointer to the NC_Dim* vector;
 * @param nvar number of NC_var;
 * @param var_vec_p pointer to the NC_Var* vector;
 * @return file pointer to a new NC_File.
 */
NC_File* nc_file_create(const char *name, int ndim, NC_Dim **dim_vec_p, int nvar, NC_Var **var_vec_p){

    NC_File *file = (NC_File*) calloc(1, sizeof(NC_File));
    /* check name length */
    file->name  = (char *) calloc(MAX_NAME_LENGTH, sizeof(char));

    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* assignment */
    file->procid = procid; /* process ncid */
    file->nprocs = nprocs; /* process number */
    /* name */
    if(snprintf(file->name, NC_MAX_NAME_LEN, "%s.%d-%d.nc", name, procid, nprocs) < 0 ){
        fprintf(stderr, "%s(%d): Error in creating NetCDF output file\n", __FILE__, __LINE__);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* dimension */
    file->Ndim     = ndim;
    file->dim_vec_p = (NC_Dim**) calloc(ndim, sizeof(NC_Dim*) );
    int i;
    for(i=0;i<ndim;i++) {file->dim_vec_p[i] = dim_vec_p[i];}
    /* var */
    file->Nvar     = nvar;
    file->var_vec_p = (NC_Var**) calloc(nvar, sizeof(NC_Var*) );
    for(i=0;i<nvar;i++) {file->var_vec_p[i] = var_vec_p[i];}

    return file;
}

/**
 * @brief read NC_Dim structure from NetCDF file.
 * @param ncid
 * @param dimid
 * @return
 */
NC_Dim* nc_dim_read_from_file(int ncid, int dimid){
    NC_Dim *dim = (NC_Dim*)calloc(1, sizeof(NC_Dim));
    dim->id = dimid;
    dim->name = calloc(NC_MAX_NAME_LEN, sizeof(char));
    nc_error( ncmpi_inq_dimname(ncid, dimid, dim->name) );
    MPI_Offset lenp;
    ncmpi_inq_dimlen(ncid, dimid, &lenp );
    dim->len = (int)lenp;
    return dim;
}

NC_Var* nc_var_read_from_file(int ncid, int varid, int NdimTol, NC_Dim **dim_p){
    NC_Var *var = (NC_Var*) calloc(1, sizeof(NC_Var));
    var->id = varid;
    int ndim;
    ncmpi_inq_varndims(ncid, varid, &(ndim));
    var->ndim = ndim;
    var->dim_vec_p = (NC_Dim**)calloc(ndim, sizeof(NC_Dim*) );
    ncmpi_inq_vartype(ncid, varid, &(var->type));
    var->name = calloc(NC_MAX_NAME_LEN, sizeof(char));
    ncmpi_inq_varname(ncid, varid, var->name);
    int n,m, dimid[ndim];
    ncmpi_inq_vardimid(ncid, varid, dimid);
    for(n=0;n<ndim;n++){
        for(m=0;m<NdimTol;m++){
            if(dimid[n] == dim_p[m]->id){ var->dim_vec_p[n] = dim_p[m]; }
        }
    }

    return var;
}
/**
 * @brief read the information from NetCDF file.
 * @param name NetCDF file name;
 * @param procid process id;
 * @param nprocs number of process;
 * @return file pointer to the NC_File structure;
 */
NC_File* nc_file_read_from_file(const char *name, int procid, int nprocs){
    NC_File *file = (NC_File*) calloc(1, sizeof(NC_File));
    /* open file */
    file->name = calloc(NC_MAX_NAME_LEN, sizeof(char));
    strcpy(file->name, name);
    nc_error( ncmpi_open(MPI_COMM_WORLD, file->name, NC_NOWRITE, MPI_INFO_NULL, &(file->ncid)) );
    /* read dimensions */
    int Ndim,n;
    nc_error( ncmpi_inq_ndims(file->ncid, &Ndim) );
    file->Ndim = Ndim;
    file->dim_vec_p = (NC_Dim**) calloc(Ndim, sizeof(NC_Dim*));
    for(n=0;n<Ndim;n++){
        file->dim_vec_p[n] = nc_dim_read_from_file(file->ncid, n);
    }
    /* read variables */
    int Nvar;
    nc_error( ncmpi_inq_nvars(file->ncid, &Nvar) );
    file->Nvar = Nvar;
    file->var_vec_p = (NC_Var**) calloc(Nvar, sizeof(NC_Var*));
    for(n=0;n<Nvar;n++){
        file->var_vec_p[n] = nc_var_read_from_file(file->ncid, n, Ndim, file->dim_vec_p);
    }
    return file;
}

/**
 * @brief print the information of nc_file.
 * @param file pointer to a nc_file structure;
 */
void nc_file_print(NC_File *file){
    printf(HEADEND "NetCDF files:\n");
    printf(HEADLINE "  name: %s\n", file->name);
    printf(HEADLINE "  procid: %d\n", file->procid);
    printf(HEADLINE "  nprocs: %d\n", file->nprocs);
    printf(HEADLINE "  Nvar: %d\n", file->Nvar);
    printf(HEADLINE "  Ndim: %d\n", file->Ndim);

    printf(HEADEND "Including variables:\n");
    int i;
    for(i=0;i<file->Nvar;i++){
        printf(HEADLINE "  name: %s\n", file->var_vec_p[i]->name);
    }
    printf(HEADEND "Including dimensions:\n");
    for(i=0;i<file->Ndim;i++){
        printf(HEADLINE "  name: %s\n", file->dim_vec_p[i]->name);
    }
}
/**
 * @brief free the memory of nc_file and associated nc_dim and nc_var.
 * @param file pointer to nc_file;
 */
void nc_file_free(NC_File *file){
    int i;
    free(file->name);
    for (i=0;i<file->Ndim;i++){ nc_dim_free(file->dim_vec_p[i]); }
    for (i=0;i<file->Nvar;i++){ nc_var_free(file->var_vec_p[i]); }
    free(file);
    return;
}

/**
 * @brief create the NetCDF file and define its dimensions and variables
 * @param file nc_file pointer
 */
void nc_file_define(NC_File *file){

    const int Ndim = file->Ndim;
    const int Nvar = file->Nvar;
    /* create output file */
    nc_error(ncmpi_create(MPI_COMM_SELF, file->name, NC_CLOBBER|NC_64BIT_OFFSET,
                          MPI_INFO_NULL, &(file->ncid)));
    /* define dimensions */
    int i,j;
    for (i=0;i<Ndim;i++){
        NC_Dim *dim = file->dim_vec_p[i]; /* pointer of dim */
        nc_error( ncmpi_def_dim(file->ncid, dim->name, dim->len, &dim->id) );
    }
    /* define variables */
    for (i=0;i<Nvar; i++){
        NC_Var *var = file->var_vec_p[i]; /* pointer of var */
        int ndim = var->ndim;
        int dimids[ndim];
        for(j=0;j<ndim;j++){ dimids[j] = var->dim_vec_p[j]->id; }
        nc_error(ncmpi_def_var(file->ncid, var->name, var->type, ndim, dimids, &(var->id) ));
    }
    /* finish definition */
    nc_error( ncmpi_enddef(file->ncid) );
}

/**
 * @brief close the NetCDF file
 * @param file nc_file pointer
 */
void nc_file_close(NC_File *file){
    nc_error( ncmpi_close(file->ncid) );

}