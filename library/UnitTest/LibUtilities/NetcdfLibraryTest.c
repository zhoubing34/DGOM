#include "LibUtiltiesTest.h"
#include "LibUtilities/NetcdfLibrary.h"

int main(int argc, char **argv){

    int procid, nprocs;

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int npoint = 10, nele = 2;
    NcDim *np = DefineNcDim("np", npoint);
    NcDim *ne = DefineNcDim("ne", nele);
    NcDim *t  = DefineNcDim("t", 0);

    NcDim **dimarray;
    NcVar **vararray;

    int ndim = 2;
    dimarray = (NcDim**) calloc(ndim, sizeof(NcDim*));
    dimarray[0] = np;
    dimarray[1] = ne;
    NcVar *x  = DefineNcVar("x", ndim, dimarray, "double");
    NcVar *y  = DefineNcVar("y", ndim, dimarray, "double");
    free(dimarray);

    ndim = 1;
    dimarray = (NcDim**) calloc(ndim, sizeof(NcDim*));
    dimarray[0] = t;
    NcVar *time = DefineNcVar("time", ndim, dimarray, "double");

    ndim = 3;
    dimarray = (NcDim**) calloc(ndim, sizeof(NcDim*));
    dimarray[0] = np;
    dimarray[1] = ne;
    dimarray[2] = t;
    int nvar = 3;
    vararray = (NcVar**) calloc(nvar, sizeof(NcVar*));
    vararray[0] = x;
    vararray[1] = y;
    vararray[2] = time;
    NcFile *file = DefineNcFile("Test", procid, nprocs, ndim, dimarray, nvar, vararray);

    CreatNcFile(file);
    CloseNcFile(file);

    MPI_Finalize();

    return 0;
}