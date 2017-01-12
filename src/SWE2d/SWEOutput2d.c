#include "swe_dirver2d.h"

/**
 * @brief
 *
 * @details
 *
 * @param[PhysDomain2d*] phys   pointer to phys structure
 * @param[SWE_Solver2d]  solver SWE solver structure
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * file     | NcFile*  | pointer to NetCDF file structure
 *
 *
 */
NcFile* SWE_SetNcOutput2d(PhysDomain2d *phys, SWE_Solver2d *solver){
    MultiReg2d   *mesh  = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    /* define dimensions */
    NcDim *ne = NcDim_create("ne", mesh->K);
    NcDim *np = NcDim_create("np", shape->Np);
    NcDim *t  = NcDim_create("t", 0);

    /* define variables */
    NcDim **dimarray;
    NcVar **vararray;

    int ndim    = 2;
    int maxDimNum = 3;
    dimarray    = (NcDim**) calloc(maxDimNum, sizeof(NcDim*));
    dimarray[0] = ne;
    dimarray[1] = np; /* the inner loop dimension comes the last */
    NcVar *x   = NcVar_create("x", ndim, dimarray, "double");
    NcVar *y   = NcVar_create("y", ndim, dimarray, "double");
    NcVar *bot = NcVar_create("bot", ndim, dimarray, "double");
//    free(dimarray);

    ndim        = 1;
//    dimarray    = (NcDim**) calloc(ndim, sizeof(NcDim*));
    dimarray[0] = t;
    NcVar *time = NcVar_create("time", ndim, dimarray, "double");
//    free(dimarray);

    ndim        = 3;
//    dimarray    = (NcDim**) calloc(ndim, sizeof(NcDim*));
    dimarray[0] = t;    /* inner loop */
    dimarray[1] = ne;
    dimarray[2] = np;
    NcVar *h   = NcVar_create("h", ndim, dimarray, "double");
    NcVar *qx  = NcVar_create("qx", ndim, dimarray, "double");
    NcVar *qy  = NcVar_create("qy", ndim, dimarray, "double");
//    free(dimarray);

    /* define files */
    ndim        = 3;
//    dimarray    = (NcDim**) calloc(ndim, sizeof(NcDim*));
    dimarray[0] = np;
    dimarray[1] = ne;
    dimarray[2] = t;
    int nvar = 7;
    vararray = (NcVar**) calloc(nvar, sizeof(NcVar*));
    vararray[0] = x;    vararray[1] = y;
    vararray[2] = time; vararray[3] = h;
    vararray[4] = qx;   vararray[5] = qy;
    vararray[6] = bot;
    NcFile *file = NcFile_create("SWE2d.", mesh->procid, mesh->nprocs, ndim, dimarray, nvar, vararray);
    free(dimarray);
    free(vararray);

    /* create files */
    NcFile_init(file);

    /* set coordinate */
    int ret;
    ret = ncmpi_put_var_double_all(file->id, file->vararray[0]->id, mesh->x[0]); NC_ERROR;
    ret = ncmpi_put_var_double_all(file->id, file->vararray[1]->id, mesh->y[0]); NC_ERROR;

    /* set bottom topography */
    ret = ncmpi_put_var_double_all(file->id, file->vararray[6]->id, solver->bot[0]); NC_ERROR;
    return file;
}

/**
 * @brief
 * Save the results into NetCDF files.
 *
 * @details
 * The results `phys->f_Q` first transfers into a float variable `var`,
 * then var is stored with `ncmpi_put_vara_float_all` function.
 *
 * @param[NcFile]       file    NetCDF file structure
 * @param[PhysDomain2d] phys    pointer to phys structure
 * @param[int]          outStep # of output steps
 * @param[double]       time    output time
 *
 * @attention
 * All variables is transfered into `float` type and then saved, while in the
 * calculation process the variables can be either `double` or `float`.
 *
 */
void SWE_StoreVar2d(NcFile *file, PhysDomain2d *phys, int outStep, double time){

    int ret;
    MPI_Offset start_v[3], count_v[3];
    MPI_Offset start_t, count_t;
    MultiReg2d *mesh  = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    /* put time */
    start_t = outStep; // start index
    count_t = 1;       // length
    ret = ncmpi_put_vara_double_all(file->id, file->vararray[2]->id, &start_t, &count_t, &time);
    NC_ERROR;

    /* put var h */
    start_v[0] = outStep;   start_v[1] = 0;         start_v[2] = 0;
    count_v[0] = 1;         count_v[1] = mesh->K;   count_v[2] = shape->Np;

    float *var = (float *) malloc(mesh->K*shape->Np*sizeof(float));
    /* put variable h */
    int k,i,ind,sk=0, Nfields=phys->Nfields;
    for (k=0;k<mesh->K;k++){
        for (i=0;i<shape->Np;i++){
            ind = (k*shape->Np + i)*Nfields;
            // printf("k = %d, i = %d, ind = %d\n", k, i, ind);
            var[sk++] = (float) phys->f_Q[ind];
        }
    }
    ret = ncmpi_put_vara_float_all(file->id, file->vararray[3]->id, start_v, count_v, var);
    NC_ERROR;

    /* put variable qx */
    sk = 0;
    for (k=0;k<mesh->K;k++){
        for (i=0;i<shape->Np;i++){
            ind = (k*shape->Np + i)*Nfields+1;
            var[sk++] = (float) phys->f_Q[ind];
        }
    }
    ret = ncmpi_put_vara_float_all(file->id, file->vararray[4]->id, start_v, count_v, var);
    NC_ERROR;

    /* put variable qy */
    sk = 0;
    for (k=0;k<mesh->K;k++){
        for (i=0;i<shape->Np;i++){
            ind = (k*shape->Np + i)*Nfields+2;
            var[sk++] = (float) phys->f_Q[ind];
        }
    }
    ret = ncmpi_put_vara_float_all(file->id, file->vararray[5]->id, start_v, count_v, var);
    NC_ERROR;

}

