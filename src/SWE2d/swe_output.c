#include <MultiRegions/Mesh/dg_mesh.h>
#include <PhysField/pf_phys.h>
#include "swe_output.h"

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
nc_file* swe_output(swe_solver *solver){

    physField *phys = solver->phys;

    const int K = phys->grid->K;
    const int Np = phys->cell->Np;
    const int procid = phys->mesh->procid;
    const int nprocs = phys->mesh->nprocs;

    /* define dimensions */
    nc_dim *ne = nc_dim_create("ne", K);
    nc_dim *np = nc_dim_create("np", Np);
    nc_dim *t  = nc_dim_create("t", 0);

    /* define variables */
    nc_dim **dimarray;
    nc_var **vararray;

    int ndim    = 2;
    size_t maxDimNum = 3;
    dimarray    = (nc_dim**) calloc(maxDimNum, sizeof(nc_dim*));
    dimarray[0] = ne;
    dimarray[1] = np; /* the inner loop dimension comes the last */
    nc_var *x   = nc_var_create("x", ndim, dimarray, NC_DOUBLE);
    nc_var *y   = nc_var_create("y", ndim, dimarray, NC_DOUBLE);
    nc_var *bot = nc_var_create("bot", ndim, dimarray, NC_DOUBLE);
//    free(dim_vec_p);

    ndim        = 1;
//    dim_vec_p    = (NcDim**) calloc(ndim, sizeof(NcDim*));
    dimarray[0] = t;
    nc_var *time = nc_var_create("time", ndim, dimarray, NC_DOUBLE);
//    free(dim_vec_p);

    ndim        = 3;
//    dim_vec_p    = (NcDim**) calloc(ndim, sizeof(NcDim*));
    dimarray[0] = t;    /* inner loop */
    dimarray[1] = ne;
    dimarray[2] = np;
    nc_var *h   = nc_var_create("h", ndim, dimarray, NC_DOUBLE);
    nc_var *qx  = nc_var_create("qx", ndim, dimarray, NC_DOUBLE);
    nc_var *qy  = nc_var_create("qy", ndim, dimarray, NC_DOUBLE);
//    free(dim_vec_p);

    /* define files */
    ndim        = 3;
//    dim_vec_p    = (NcDim**) calloc(ndim, sizeof(NcDim*));
    dimarray[0] = np;
    dimarray[1] = ne;
    dimarray[2] = t;
    int nvar = 7;
    vararray = (nc_var**) calloc((size_t)nvar, sizeof(nc_var*));
    vararray[0] = x;    vararray[1] = y;
    vararray[2] = time; vararray[3] = h;
    vararray[4] = qx;   vararray[5] = qy;
    vararray[6] = bot;
    nc_file *file = nc_file_create("swe2d.", procid, nprocs, ndim, dimarray, nvar, vararray);
    free(dimarray);
    free(vararray);

    /* create files */
    nc_file_init(file);

    /* set coordinate */
    ncmpi_put_var_double_all(file->id, file->var_vec_p[0]->id, phys->region->x[0]);
    ncmpi_put_var_double_all(file->id, file->var_vec_p[1]->id, phys->region->y[0]);

    /* set bottom topography */
    ncmpi_put_var_double_all(file->id, file->var_vec_p[6]->id, solver->bot);
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
 * @param[in] solver SWE solver structure
 * @param[in] tstep # of output steps
 * @param[in] time output time
 *
 * @attention
 * All variables is transfered into `float` type and then saved, while in the
 * calculation process the variables can be either `double` or `float`.
 *
 */
void swe_save_var(swe_solver *solver, int tstep, double time){

    physField *phys = solver->phys;

    const int K = phys->grid->K;
    const int Np = phys->cell->Np;
    nc_file *file = solver->outfile;

    MPI_Offset start_v[3], count_v[3];
    MPI_Offset start_t, count_t;

    /* put time */
    start_t = tstep; // start index
    count_t = 1;       // length
    ncmpi_put_vara_double_all(file->id, file->var_vec_p[2]->id, &start_t, &count_t, &time);

    /* put var h */
    start_v[0] = tstep;
    start_v[1] = 0;
    start_v[2] = 0;
    count_v[0] = 1;
    count_v[1] = K;
    count_v[2] = Np;

    float *var = (float *) malloc(K*Np*sizeof(float));
    /* put variable h */
    int k,i,ind,sk=0, Nfields = phys->Nfield;
    for (k=0;k< K;k++){
        for (i=0;i<Np;i++){
            ind = (k*Np + i)*Nfields;
            var[sk++] = (float) phys->f_Q[ind];
        }
    }
    ncmpi_put_vara_float_all(file->id, file->var_vec_p[3]->id, start_v, count_v, var);

    /* put variable qx */
    sk = 0;
    for (k=0;k<K;k++){
        for (i=0;i<Np;i++){
            ind = (k*Np + i)*Nfields+1;
            var[sk++] = (float) phys->f_Q[ind];
        }
    }
    ncmpi_put_vara_float_all(file->id, file->var_vec_p[4]->id, start_v, count_v, var);

    /* put variable qy */
    sk = 0;
    for (k=0;k<K;k++){
        for (i=0;i<Np;i++){
            ind = (k*Np + i)*Nfields+2;
            var[sk++] = (float) phys->f_Q[ind];
        }
    }
    ncmpi_put_vara_float_all(file->id, file->var_vec_p[5]->id, start_v, count_v, var);

}

