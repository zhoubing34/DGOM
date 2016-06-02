#include <pnetcdf.h>

static void handle_error(int status, int lineno)
{
    fprintf(stderr, "Error at line %d: %s\n", lineno, ncmpi_strerror(status));
    MPI_Abort(MPI_COMM_WORLD, 1);
}


typedef struct ncstru {
    int ncfile;     // file handle
    int * varid;    // variable id
    char * filename; // filename
}Ncfile;


Ncfile* SetupOutput(Mesh *, char* );

void PutVar(Ncfile * , int , double , Mesh * );