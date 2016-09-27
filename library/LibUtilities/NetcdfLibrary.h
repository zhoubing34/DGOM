#ifndef NCLIBRARY_H
#define NCLIBRARY_H

#include "pnetcdf.h"

/* NetcdfLibrary */
#define NCNAMELEN 1024

typedef struct {
    char *name; /* name of dimension */
    int  len;   /* length of dimension */
    int  id;    /* dimension id in NetCDF file */
}NcDim;

typedef struct {
    char  *name;        /* name of variables */
    int   ndim;
    NcDim **dimarray;   /* array of dimension pointers */
    int   type;
    int   id;           /* variable id in NetCDF file */
}NcVar;

typedef struct {
    char  *name;        /* name of files */
    int   procid;
    int   nprocs;
    int   ndim, nvar;
    NcDim **dimarray;   /* array of dimension pointers */
    NcVar **vararray;   /* array of variable pointers */
    int   id;           /* file id of NetCDF file */
}NcFile;

/* handing error marcos */
#define NC_ERROR  \
do{ \
    if(ret != NC_NOERR){ \
        handle_error(ret, __LINE__); \
    } \
}while(0)

/* public function */
void handle_error(int status, int lineno);
NcDim*  NcDim_create(char *namestr, int len);
NcVar*  NcVar_create(char *namestr, int ndim, NcDim **dimArray, char *typestr);
NcFile* NcFile_create(char *namestr, int procid, int nprocs,
                      int ndim, NcDim **dimArray,
                      int nvar, NcVar **varArray);
void NcFile_init(NcFile *file);
void NcFile_close(NcFile *file);

void NcFile_free(NcFile *file);

void NcDim_print(NcDim *dim);
void NcVar_print(NcVar *var);
void NcFile_print(NcFile *file);

#endif