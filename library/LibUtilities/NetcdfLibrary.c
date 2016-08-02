#include "NetcdfLibrary.h"
#include "LibUtilities.h"

NcDim* DefineNcDim(char *namestr, int len){
    NcDim *dim = (NcDim*) calloc(1, sizeof(NcDim));
    dim->name  = (char *) malloc(NCNAMELEN*sizeof(char));

    /* check name length */
    size_t namelen = strlen(namestr);
    if (namelen > NCNAMELEN){
        printf("The name of dimension %s is too long!\n", namestr);
        exit(-1);
    }

    /* assignment */
    strcpy(dim->name, namestr);
    if (len>0) {
        dim->len = len;
    }else if (len==0){
        dim->len = NC_UNLIMITED;
    }

    return dim;
}

void FreeNcDim(NcDim *dim){
    free(dim->name);
    free(dim);
}

void PrintNcDim(NcDim *dim){
    printf("Netcdf Dimension\n");
    printf("name: %s\n", dim->name);
    printf("length: %d\n", dim->len);
}

NcVar* DefineNcVar(char *namestr, int ndim, NcDim **dimArray, char *typestr){
    NcVar *var = (NcVar*) calloc(1, sizeof(NcVar));
    var->name  = (char *) malloc(NCNAMELEN*sizeof(char));

    /* check name length */
    size_t namelen = strlen(namestr);
    if (namelen > NCNAMELEN){
        printf("The name of variable %s is too long!\n", namestr);
        exit(-1);
    }

    /* assignment */
    strcpy(var->name, namestr);
    var->ndim     = ndim; /* number of dimensions */
    var->dimarray = (NcDim**) calloc(ndim, sizeof(NcDim*) );
    int i;
    for(i=0;i<ndim;i++)
        var->dimarray[i] = dimArray[i];

    if ( !(memcmp(typestr, "int", 3)) ){
        var->type = NC_INT;
    }else if( !(memcmp(typestr, "float", 5)) ){
        var->type = NC_FLOAT;
    }else if( !(memcmp(typestr, "double", 6)) ){
        var->type = NC_DOUBLE;
    }else if( !(memcmp(typestr, "short", 5)) ){
        var->type = NC_SHORT;
    }else if( !(memcmp(typestr, "string", 6)) ){
        var->type = NC_STRING;
    }else{
        printf("Wrong type of NetCDF variable: %s!\n", typestr);
        exit(-1);
    }
    return var;
}

void PrintNcVar(NcVar *var){
    printf("NetCDF Variable\n");
    printf("name: %s\n", var->name);
    printf("ndim: %d\n", var->ndim);
    printf("type: %d\n", var->type);
    printf("Dimensions contained\n");

    int i;
    for(i=0;i<var->ndim;i++){
        PrintNcDim(var->dimarray[i]);
    }
}

void FreeNcVar(NcVar *var){
    free(var->name);
    free(var->dimarray);
    free(var);
}

NcFile* DefineNcFile(char *namestr, int procid, int nprocs,
                     int ndim, NcDim **dimArray,
                     int nvar, NcVar **varArray){

    NcFile *file = (NcFile*) calloc(1, sizeof(NcFile));
    file->name  = (char *) malloc(NCNAMELEN*sizeof(char));

    /* check name length */
    size_t namelen = strlen(namestr);
    if (namelen > NCNAMELEN){
        printf("The name of variable %s is too long!\n", namestr);
        exit(-1);
    }

    /* assignment */
    file->procid = procid; /* process id */
    file->nprocs = nprocs; /* process number */
    /* name */
    int ret = snprintf(file->name, NCNAMELEN, "%s%d-%d.nc", namestr, procid, nprocs);

    /* dimension */
    file->ndim     = ndim;
    file->dimarray = (NcDim**) calloc(ndim, sizeof(NcDim*) );
    int i;
    for(i=0;i<ndim;i++)
        file->dimarray[i] = dimArray[i];
    /* var */
    file->nvar     = nvar;
    file->vararray = (NcVar**) calloc(nvar, sizeof(NcVar*) );
    for(i=0;i<nvar;i++)
        file->vararray[i] = varArray[i];

    return file;
}

void PrintNcFile(NcFile *file){
    printf("NetCDF Files\n");
    printf("name: %s\n", file->name);
    printf("procid: %d\n", file->procid);
    printf("nprocs: %d\n", file->nprocs);
    printf("nvar: %d\n", file->nvar);
    printf("ndim: %d\n", file->ndim);

    int i;
    for(i=0;i<file->nvar;i++){
        PrintNcVar(file->vararray[i]);
    }
    for(i=0;i<file->ndim;i++){
        PrintNcDim(file->dimarray[i]);
    }
}

void FreeNcFile(NcFile *file){
    int i;
    free(file->name);
    for (i=0;i<file->ndim;i++)
        FreeNcDim(file->dimarray[i]);

    for (i=0;i<file->nvar;i++)
        FreeNcVar(file->vararray[i]);

    free(file);
}

void CreatNcFile(NcFile *file){

    int ret;
    /* create output file */
    ret = ncmpi_create(MPI_COMM_SELF, file->name, NC_CLOBBER|NC_64BIT_OFFSET, MPI_INFO_NULL, &file->id);
    NC_ERROR;

    /* define dimensions */
    int i;
    NcDim *dim;
    for (i=0;i<file->ndim;i++){
        dim = file->dimarray[i]; /* pointer of dim */

        ret = ncmpi_def_dim(file->id, dim->name, dim->len, &dim->id);
        NC_ERROR;
    }

    /* define variables */
    int *dimids, j;
    NcVar *var;
    for (i=0; i<file->nvar; i++){
        var    = file->vararray[i]; /* pointer of var */
        dimids = (int *)calloc(var->ndim, sizeof(int));
        for (j=0; j<var->ndim; j++)
            dimids[j] = var->dimarray[j]->id;

        ret = ncmpi_def_var(file->id, var->name, var->type, var->ndim, dimids, &var->id);
        NC_ERROR;
        free(dimids);
    }

    /* finish definition */
    ret = ncmpi_enddef(file->id); NC_ERROR;
}

void CloseNcFile(NcFile *file){
    int ret;
    ret = ncmpi_close(file->id);
    NC_ERROR;
}

/**
 * @brief pnetcdf error handle function
 * @details
 * @param[in] status error number
 * @param[in] lineno line number
 */
void handle_error(int status, int lineno){
    fprintf(stderr, "Error at line %d: %s\n", lineno, ncmpi_strerror(status));
    MPI_Abort(MPI_COMM_WORLD, 1);
}