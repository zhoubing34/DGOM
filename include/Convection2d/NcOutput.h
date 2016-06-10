/** \file   NcOutput.c
 *  \brief
 *  pnetcdf output function declarations
 */

#ifndef NCOUTPUT_H
#define NCOUTPUT_H

#include <pnetcdf.h>

/**
 * @brief pnetcdf error handle function
 * @details
 * @param[in] status error number
 * @param[in] lineno line number
 */
static void handle_error(int status, int lineno)
{
    fprintf(stderr, "Error at line %d: %s\n", lineno, ncmpi_strerror(status));
    MPI_Abort(MPI_COMM_WORLD, 1);
}

/* check if it returns error */
#define NC_ERROR if(ret != NC_NOERR) handle_error(ret, __LINE__);

typedef struct Ncfile {
    /** nc file handle */
    int ncfile;
    /** array of variable id */
    int * varid;
    /** nc file name */
    char * filename;
}Ncfile;

/* NcOutput.c */
Ncfile* SetupOutput(Mesh *, char* );
void PutVar(Ncfile * , int , double , Mesh * );

#endif