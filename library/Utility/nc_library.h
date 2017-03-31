/**
 * @file
 * NetCDF相关结构体定义
 * @details
 * 1. `nc_dim_create` generates a pointer to a new nc_dim;
 * 2. `nc_var_create` generates a pointer to a new nc_var;
 * 3. `nc_file_create` generates a pointer to a new nc_file;
 * 4. `nc_file_init` creates the NetCDF file and defines its dimensions and variables;
 * 5. `nc_file_free` free the memory of nc_file and its including nc_dim and nc_var;
 * 6. `nc_file_close` close the NetCDF file associated with the nc_file.
 * 7. `nc_file_print` print all the information of nc_file.
 */
#ifndef NCLIBRARY_H
#define NCLIBRARY_H

#include "pnetcdf.h"

/** @brief structure for NetCDF dimensions */
typedef struct {
    char *name; ///< name of dimension
    int  len; ///< length of dimension
    int  id; ///< dimension id in NetCDF file
}NC_Dim;
/** @brief structure for NetCDF variable */
typedef struct {
    char *name; ///< name of variables
    int ndim; ///< number of dimension in the variable
    NC_Dim **dim_vec_p; ///< pointer to the array of nc_dim
    int type; ///< variable type
    int id; ///< variable id in NetCDF file
}NC_Var;
/** @brief structure for NetCDF file */
typedef struct {
    char *name; ///< name of files;
    int procid; ///< process id;
    int nprocs; ///< number of processes;
    int Ndim; ///< number of dimensions;
    int Nvar; ///< number of variable;
    NC_Dim **dim_vec_p; ///< pointer to the array of nc_dim;
    NC_Var **var_vec_p; ///< pointer to the array of nc_var;
    int ncid; ///< file id of NetCDF file;
}NC_File;

/* public function */
NC_Dim* nc_dim_create(const char *name, int len);
NC_Var* nc_var_create(const char *name, int ndim, NC_Dim **dim_vec_p, int type);
NC_File* nc_file_create(const char *name, int ndim, NC_Dim **dim_vec_p, int nvar, NC_Var **var_vec_p);
void nc_file_define(NC_File *file);
NC_File* nc_file_read_from_file(const char *name, int procid, int nprocs);
void nc_file_close(NC_File *file);
void nc_file_free(NC_File *file);
void nc_file_print(NC_File *file);
void nc_var_print(NC_Var *var);

/** marco for handing error of parallel-netcdf */
#define nc_error(t)                                         \
do{                                                         \
    if(t != NC_NOERR){                                      \
        fprintf(stderr, "%s: error at line %d: %s\n",       \
                __FUNCTION__, __LINE__, ncmpi_strerror(t)); \
        MPI_Abort(MPI_COMM_WORLD, 1);                       \
    }                                                       \
}while(0)

#endif