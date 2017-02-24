/**
 * @file
 * utility functions
 * @details
 * 1. `str2int` transfer string to integer;
 * 2. `str2double` transfer string to double;
 * 3. utility function in Utility library.
 */

#ifndef LIBUTILITIES_H
#define LIBUTILITIES_H

#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "pnetcdf.h"

/* float type */
#ifdef DOUBLE_PRECISION
typedef double dg_real; ///< user specific float type
#define MPI_TYPE MPI_DOUBLE ///< variable type for MPI subroutines
#define NC_TYPE NC_DOUBLE ///< variable type for netCDF subroutines
#else
typedef float dg_real; ///< user specific float type
#define MPI_TYPE MPI_FLOAT ///< variable type for MPI subroutines
#define NC_TYPE NC_FLOAT ///< variable type for netCDF subroutines
#endif

/* max character length */
#define MAX_NAME_LENGTH 1024
#define EPS 100*DBL_EPSILON

#include "mat_utils.h"

/* file I/O */
/* open file and print error message and exit if fails */
#define dg_fopen(fp, filename, errmessage) \
if( (fp = fopen(filename, "r")) == NULL ){ \
    fprintf(stderr, "%s (%s): %d\n%s: %s.\n", \
            __FUNCTION__, __FILE__,__LINE__,errmessage,filename);\
    exit(-1); \
}

/* string to int */
void str2int(char *str, int *N, char* errmessage);
void str2double(char *str, double *scal, char* errmessage);

#endif //LIBUTILITIES_H