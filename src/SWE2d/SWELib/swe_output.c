#include "swe_lib.h"

/**
 * @brief
 */
void swe_output(){

    extern SWE_Solver solver;
    dg_phys *phys = solver.phys;

    const int K = dg_grid_K( dg_phys_grid(phys) );
    const int Np = dg_cell_Np( dg_phys_cell(phys) );
    const int procid = dg_grid_procid( dg_phys_grid(phys) );
    const int nprocs = dg_grid_nprocs( dg_phys_grid(phys) );

    /* define dimensions */
    NC_Dim *ne = nc_dim_create("ne", K);
    NC_Dim *np = nc_dim_create("np", Np);
    NC_Dim *t  = nc_dim_create("t", 0);

    /* define variables */
    NC_Dim **dimarray;
    NC_Var **vararray;

    int ndim    = 2;
    size_t maxDimNum = 3;
    dimarray    = (NC_Dim**) calloc(maxDimNum, sizeof(NC_Dim*));
    dimarray[0] = ne;
    dimarray[1] = np; /* the inner loop dimension comes the last */
    NC_Var *x   = nc_var_create("x", ndim, dimarray, NC_FLOAT);
    NC_Var *y   = nc_var_create("y", ndim, dimarray, NC_FLOAT);
    NC_Var *bot = nc_var_create("bot", ndim, dimarray, NC_FLOAT);

    ndim        = 1;
    dimarray[0] = t;
    NC_Var *time = nc_var_create("time", ndim, dimarray, NC_FLOAT);

    ndim        = 3;
    dimarray[0] = t;
    dimarray[1] = ne;
    dimarray[2] = np;
    NC_Var *h   = nc_var_create("h", ndim, dimarray, NC_FLOAT);
    NC_Var *qx  = nc_var_create("qx", ndim, dimarray, NC_FLOAT);
    NC_Var *qy  = nc_var_create("qy", ndim, dimarray, NC_FLOAT);

    /* define files */
    ndim        = 3;
    dimarray[0] = np;
    dimarray[1] = ne;
    dimarray[2] = t;
    int nvar = 7;
    vararray = (NC_Var**) calloc((size_t)nvar, sizeof(NC_Var*));
    vararray[0] = x;    vararray[1] = y;
    vararray[2] = time; vararray[3] = h;
    vararray[4] = qx;   vararray[5] = qy;
    vararray[6] = bot;
    NC_File *file = nc_file_create("swe2d.", procid, nprocs, ndim, dimarray, nvar, vararray);
    free(dimarray);
    free(vararray);

    /* create files */
    nc_file_define(file);

    /* coordinate x */
    const int Nfields = dg_phys_Nfield(phys);
    int k,n;
    float *var = (float *) malloc(K*Np*sizeof(float));

    for(k=0;k<K;k++){
        for(n=0;n<Np;n++){
            var[k*Np+n] = (float) dg_region_x(dg_phys_region(phys))[k][n];
        }
    }
    ncmpi_put_var_float_all(file->ncid, file->var_vec_p[0]->id, var);

    /* coordinate y */
    for(k=0;k<K;k++){
        for(n=0;n<Np;n++){
            var[k*Np+n] = (float) dg_region_y(dg_phys_region(phys))[k][n];
        }
    }
    ncmpi_put_var_float_all(file->ncid, file->var_vec_p[1]->id, var);

    /* bottom elevation */
    for (k=0;k<K;k++){
        for (n=0;n<Np;n++){
            int ind = (k*Np + n)*Nfields+3;
            var[k*Np+n] = (float) dg_phys_f_Q(phys)[ind];
        }
    }
    ncmpi_put_var_float_all(file->ncid, file->var_vec_p[6]->id, var);

    /* free the memory */
    free(var);
    /* assignment */
    solver.outfile = file;
    return;
}

/**
 * @brief
 * Save the results into NetCDF files.
 *
 * @details
 * The results `phys->f_Q` first transfers into a float variable `var`,
 * then var is stored with `ncmpi_put_vara_float_all` function.
 *
 * @param[in] tstep # of output steps
 * @param[in] time output time
 *
 * @attention
 * All variables is transfered into `float` type and then saved, while in the
 * calculation process the variables can be either `double` or `float`.
 */
void swe_save_var(int tstep, double time){

    extern SWE_Solver solver;
    dg_phys *phys = solver.phys;

    const int K = dg_grid_K( dg_phys_grid(phys) );
    const int Np = dg_cell_Np( dg_phys_cell(phys) );
    const int Nfields = dg_phys_Nfield(phys);

    dg_real *f_Q = dg_phys_f_Q(phys);
    NC_File *file = solver.outfile;

    MPI_Offset start_v[3] = {tstep, 0, 0};
    MPI_Offset count_v[3] = {1, K, Np};
    MPI_Offset start_t = tstep;
    MPI_Offset count_t = 1;

    /* put time */
    ncmpi_put_vara_double_all(file->ncid, file->var_vec_p[2]->id, &start_t, &count_t, &time);

    /* put var h */
    float *var = (float *) malloc(K*Np*sizeof(float));
    /* put variable h */
    int k,i,ind,sk=0;
    for (k=0;k<K;k++){
        for (i=0;i<Np;i++){
            ind = (k*Np + i)*Nfields;
            var[sk++] = (float) f_Q[ind];
        }
    }
    ncmpi_put_vara_float_all(file->ncid, file->var_vec_p[3]->id, start_v, count_v, var);

    /* put variable qx */
    sk = 0;
    for (k=0;k<K;k++){
        for (i=0;i<Np;i++){
            ind = (k*Np + i)*Nfields+1;
            var[sk++] = (float) f_Q[ind];
        }
    }
    ncmpi_put_vara_float_all(file->ncid, file->var_vec_p[4]->id, start_v, count_v, var);

    /* put variable qy */
    sk = 0;
    for (k=0;k<K;k++){
        for (i=0;i<Np;i++){
            ind = (k*Np + i)*Nfields+2;
            var[sk++] = (float) f_Q[ind];
        }
    }
    ncmpi_put_vara_float_all(file->ncid, file->var_vec_p[5]->id, start_v, count_v, var);

    free(var);
    return;
}

