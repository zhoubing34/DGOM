#include "SWEDriver2d.h"

NcFile* SWEOutput(PhysDomain2d *phys, SWESolver *solver){
    MultiReg2d   *mesh  = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    /* define dimensions */
    NcDim *ne = DefineNcDim("ne", mesh->K);
    NcDim *np = DefineNcDim("np", shape->Np);
    NcDim *t  = DefineNcDim("t",  0);

    /* define variables */
    NcDim **dimarray;
    NcVar **vararray;

    int ndim = 2;
    dimarray = (NcDim**) calloc(ndim, sizeof(NcDim*));
    dimarray[0] = ne;
    dimarray[1] = np;
    NcVar *x   = DefineNcVar("x", ndim, dimarray, "double");
    NcVar *y   = DefineNcVar("y", ndim, dimarray, "double");
    NcVar *bot = DefineNcVar("bot", ndim, dimarray, "double");
    free(dimarray);

    ndim = 1;
    dimarray = (NcDim**) calloc(ndim, sizeof(NcDim*));
    dimarray[0] = t;
    NcVar *time = DefineNcVar("time", ndim, dimarray, "double");
    free(dimarray);

    ndim = 3;
    dimarray = (NcDim**) calloc(ndim, sizeof(NcDim*));
    dimarray[0] = t;    /* inner loop */
    dimarray[1] = ne;
    dimarray[2] = np;
    NcVar *h   = DefineNcVar("h",  ndim, dimarray, "double");
    NcVar *qx  = DefineNcVar("qx", ndim, dimarray, "double");
    NcVar *qy  = DefineNcVar("qy", ndim, dimarray, "double");
    free(dimarray);

    /* define files */
    ndim = 3;
    dimarray = (NcDim**) calloc(ndim, sizeof(NcDim*));
    dimarray[0] = np;
    dimarray[1] = ne;
    dimarray[2] = t;
    int nvar = 7;
    vararray = (NcVar**) calloc(nvar, sizeof(NcVar*));
    vararray[0] = x;    vararray[1] = y;
    vararray[2] = time; vararray[3] = h;
    vararray[4] = qx;   vararray[5] = qy;
    vararray[6] = bot;
    NcFile *file = DefineNcFile("SWE2D", mesh->procid, mesh->nprocs, ndim, dimarray, nvar, vararray);
    free(dimarray);
    free(vararray);

    /* create files */
    CreatNcFile(file);

    /* set coordinate */
    int ret;
    ret = ncmpi_put_var_double_all(file->id, file->vararray[0]->id, mesh->x[0]); NC_ERROR;
    ret = ncmpi_put_var_double_all(file->id, file->vararray[1]->id, mesh->y[0]); NC_ERROR;

    /* set bottom topography */
    ret = ncmpi_put_var_double_all(file->id, file->vararray[6]->id, solver->bot[0]); NC_ERROR;
    return file;
}


void StoreVar(NcFile *file, PhysDomain2d *phys, int outStep, double time){

    int ret;
    MPI_Offset start_v[3], count_v[3];
    MPI_Offset start_t, count_t;
    MultiReg2d   *mesh  = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    /* put time */
    start_t = outStep; // start index
    count_t = 1;       // length
    ret = ncmpi_put_vara_double_all(file->id, file->vararray[2]->id, &start_t, &count_t, &time); NC_ERROR;

    /* put var h */
    start_v[0] = outStep;
    start_v[1] = 0;
    start_v[2] = 0;
    count_v[0] = 1;
    count_v[1] = mesh->K;
    count_v[2] = shape->Np;

    real *var = (real *) malloc(mesh->K*shape->Np*sizeof(real));
    int k,i,ind,sk=0, Nfields=phys->Nfields;
    for (k=0;k<mesh->K;k++){
        for (i=0;i<shape->Np;i++){
            ind = (k*shape->Np + i)*Nfields;
            // printf("k = %d, i = %d, ind = %d\n", k, i, ind);
            var[sk++] = phys->f_Q[ind];
        }
    }
    ret = ncmpi_put_vara_float_all(file->id, file->vararray[3]->id, start_v, count_v, var); NC_ERROR;

    /* put var qx */
    sk = 0;
    for (k=0;k<mesh->K;k++){
        for (i=0;i<shape->Np;i++){
            ind = (k*shape->Np + i)*Nfields+1;
            var[sk++] = phys->f_Q[ind];
        }
    }
    ret = ncmpi_put_vara_float_all(file->id, file->vararray[4]->id, start_v, count_v, var); NC_ERROR;

    /* put var qy */
    sk = 0;
    for (k=0;k<mesh->K;k++){
        for (i=0;i<shape->Np;i++){
            ind = (k*shape->Np + i)*Nfields+2;
            var[sk++] = phys->f_Q[ind];
        }
    }
    ret = ncmpi_put_vara_float_all(file->id, file->vararray[5]->id, start_v, count_v, var); NC_ERROR;

}

