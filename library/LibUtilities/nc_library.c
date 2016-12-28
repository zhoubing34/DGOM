#include "nc_library.h"
#include "LibUtilities.h"

#define NCNAMELEN 1024

static void handle_error(int status, int lineno);

nc_dim* nc_dim_create(char *namestr, int len){
    nc_dim *dim = (nc_dim*) calloc(1, sizeof(nc_dim));

    /* allocate name */
    size_t namelen = strlen(namestr);
    dim->name  = (char *) calloc(namelen, sizeof(char));

    /* assignment */
    strcpy(dim->name, namestr);
    if (len>0) {
        dim->len = len;
    }else if (len==0){
        dim->len = NC_UNLIMITED;
    }

    return dim;
}

void nc_dim_free(nc_dim *dim){
    free(dim->name);
    free(dim);
}

void nc_dim_print(nc_dim *dim){
    printf("NetCDF Dimension\n");
    printf("  name:   %s\n", dim->name);
    printf("  length: %d\n", dim->len);
    printf("\n");
}

nc_var* nc_var_create(char *namestr, int ndim, nc_dim **dimArray, int nc_type){
    nc_var *var = (nc_var*) calloc(1, sizeof(nc_var));

    /* allocate name length */
    size_t namelen = strlen(namestr);
    var->name  = (char *) calloc(namelen, sizeof(char));

    /* assignment */
    strcpy(var->name, namestr);
    var->ndim     = ndim; /* number of dimensions */
    var->dimarray = (nc_dim**) calloc(ndim, sizeof(nc_dim*) );
    int i;
    for(i=0;i<ndim;i++)
        var->dimarray[i] = dimArray[i];

    switch (nc_type){
        case NC_INT: var->type = NC_INT; break;
        case NC_FLOAT: var->type = NC_FLOAT; break;
        case NC_DOUBLE: var->type = NC_DOUBLE; break;
        case NC_SHORT: var->type = NC_SHORT; break;
        case NC_STRING: var->type = NC_STRING; break;
        default: printf("Unknown type of NetCDF variable type: %d!\n", nc_type);
    }
    return var;
}

void nc_var_free(nc_var *var){
    free(var->name);
    free(var->dimarray);
    free(var);
}

void nc_var_print(nc_var *var){
    printf("NetCDF Variable\n");
    printf("  name: %s\n", var->name);
    printf("  ndim: %d\n", var->ndim);
    printf("  type: %d\n", var->type);
    printf("  Dimensions contained\n");

    int i;
    for(i=0;i<var->ndim;i++){
        nc_dim_print(var->dimarray[i]);
    }
    printf("\n");
}

nc_file* nc_file_create(char *namestr, int procid, int nprocs,
                       int ndim, nc_dim **dimArray,
                       int nvar, nc_var **varArray){

    nc_file *file = (nc_file*) calloc(1, sizeof(nc_file));

    /* check name length */
    size_t namelen = strlen(namestr);
    file->name  = (char *) calloc(namelen, sizeof(char));

    /* assignment */
    file->procid = procid; /* process id */
    file->nprocs = nprocs; /* process number */
    /* name */
    if(snprintf(file->name, NCNAMELEN, "%s%d-%d.nc", namestr, procid, nprocs) < 0 ){
        fprintf(stderr, "Error in writing NetCDF output file: at line %d, %s\n", __LINE__, __FILE__);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* dimension */
    file->ndim     = ndim;
    file->dimarray = (nc_dim**) calloc(ndim, sizeof(nc_dim*) );
    int i;
    for(i=0;i<ndim;i++)
        file->dimarray[i] = dimArray[i];
    /* var */
    file->nvar     = nvar;
    file->vararray = (nc_var**) calloc(nvar, sizeof(nc_var*) );
    for(i=0;i<nvar;i++)
        file->vararray[i] = varArray[i];

    return file;
}

void nc_file_print(nc_file *file){
    printf("NetCDF Files\n");
    printf("  name:   %s\n", file->name);
    printf("  procid: %d\n", file->procid);
    printf("  nprocs: %d\n", file->nprocs);
    printf("  nvar:   %d\n", file->nvar);
    printf("  ndim:   %d\n", file->ndim);

    int i;
    for(i=0;i<file->nvar;i++){
        nc_var_print(file->vararray[i]);
    }
    for(i=0;i<file->ndim;i++){
        nc_dim_print(file->dimarray[i]);
    }
}

void nc_file_free(nc_file *file){
    int i;
    free(file->name);
    for (i=0;i<file->ndim;i++)
        nc_dim_free(file->dimarray[i]);

    for (i=0;i<file->nvar;i++)
        nc_var_free(file->vararray[i]);

    free(file);
}

/* Init the NetCDF files */
void nc_file_init(nc_file *file){

    /* create output file */
    NC_ERROR(ncmpi_create(MPI_COMM_SELF, file->name, NC_CLOBBER|NC_64BIT_OFFSET, MPI_INFO_NULL, &file->id));

    /* define dimensions */
    int i;
    nc_dim *dim;
    for (i=0;i<file->ndim;i++){
        dim = file->dimarray[i]; /* pointer of dim */

        NC_ERROR( ncmpi_def_dim(file->id, dim->name, dim->len, &dim->id) );
    }

    /* define variables */
    int *dimids, j;
    nc_var *var;
    for (i=0; i<file->nvar; i++){
        var    = file->vararray[i]; /* pointer of var */
        dimids = (int *)calloc(var->ndim, sizeof(int));
        for (j=0; j<var->ndim; j++)
            dimids[j] = var->dimarray[j]->id;

        NC_ERROR(ncmpi_def_var(file->id, var->name, var->type, var->ndim, dimids, &var->id));
        free(dimids);
    }

    /* finish definition */
    NC_ERROR( ncmpi_enddef(file->id) );
}

/* close NetCDF file */
void nc_file_close(nc_file *file){

    NC_ERROR( ncmpi_close(file->id) );

}

///**
// * @brief pnetcdf error handle function
// * @details
// * @param[in] status error number
// * @param[in] lineno line number
// */
//void handle_error(int status, int lineno){
//    fprintf(stderr, "Error at line %d: %s\n", lineno, ncmpi_strerror(status));
//    MPI_Abort(MPI_COMM_WORLD, 1);
//}