#ifndef NCLIBRARY_H
#define NCLIBRARY_H

#include "pnetcdf.h"

/* NetCDF Library */

typedef struct {
    char *name; /* name of dimension */
    int  len;   /* length of dimension */
    int  id;    /* dimension id in NetCDF file */
}nc_dim;

typedef struct {
    char  *name;        /* name of variables */
    int   ndim;
    nc_dim **dimarray;   /* array of dimension pointers */
    int   type;
    int   id;           /* variable id in NetCDF file */
}nc_var;

typedef struct {
    char  *name;        /* name of files */
    int   procid;
    int   nprocs;
    int   ndim, nvar;
    nc_dim **dimarray;   /* array of dimension pointers */
    nc_var **vararray;   /* array of variable pointers */
    int   id;           /* file id of NetCDF file */
}nc_file;

/* handing error marcos */
#define NC_ERROR(t)  \
do{ \
    if(t != NC_NOERR){ \
        fprintf(stderr, "Error at line %d: %s\n", __LINE__, ncmpi_strerror(t));\
        MPI_Abort(MPI_COMM_WORLD, 1);\
    } \
}while(0)

/* public function */
nc_dim*  nc_dim_create(char *namestr, int len);
nc_var* nc_var_create(char *namestr, int ndim, nc_dim **dimArray, int nc_type);
nc_file* nc_file_create(char *namestr, int procid, int nprocs,
                       int ndim, nc_dim **dimArray,
                       int nvar, nc_var **varArray);
void nc_file_init(nc_file *file);
void nc_file_close(nc_file *file);
void nc_file_free(nc_file *file);

void nc_dim_print(nc_dim *dim);
void nc_var_print(nc_var *var);
void nc_file_print(nc_file *file);


void handle_error(int status, int lineno);
#endif