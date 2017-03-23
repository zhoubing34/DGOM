#include "conv_driver.h"

/**
 * @brief
 * Create output files in NetCDF format
 *
 * @details
 * Each process opens a file and returns its file object Ncfile
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 * @return
 * Ncfile* outfile
 * fields   | type      | description of value
 * -------- |---------- |----------------------
 * ncfile   | int       | netcdf file handle
 * varid    | int*      | variable ids
 *
 */
void conv_setoutput(){

    extern Conv_Solver solver;
    dg_phys *phys = solver.phys;
    dg_grid *grid = dg_phys_grid(phys);
    const int K = dg_grid_K(grid);
    const int Np = dg_cell_Np(dg_phys_cell(phys));
    const int procid = dg_grid_procid(grid);
    const int nprocs = dg_grid_nprocs(grid);

    /* define dimensions */
    NC_Dim *ne = nc_dim_create("ne", K);
    NC_Dim *np = nc_dim_create("np", Np);
    NC_Dim *t  = nc_dim_create("t", 0);

    const int ndim = 3;
    const int nvar = 4;
    NC_Dim *dimarray[ndim];
    NC_Var *vararray[nvar];

    /* define variables: x,y,time and h */
    dimarray[0] = ne;
    dimarray[1] = np; /* the inner loop dimension comes the last */
    NC_Var *x = nc_var_create("x", 2, dimarray, NC_DOUBLE);
    NC_Var *y = nc_var_create("y", 2, dimarray, NC_DOUBLE);

    dimarray[0] = t;
    NC_Var *time = nc_var_create("time", 1, dimarray, NC_DOUBLE);

    dimarray[0] = t;
    dimarray[1] = ne;
    dimarray[2] = np; /* inner loop */
    NC_Var *h = nc_var_create("var", 3, dimarray, NC_FLOAT);

    /* define files */
    dimarray[0] = np;
    dimarray[1] = ne;
    dimarray[2] = t;

    vararray[0] = x;
    vararray[1] = y;
    vararray[2] = time;
    vararray[3] = h;
    NC_File *file = nc_file_create("conv2d.", procid, nprocs, ndim, dimarray, nvar, vararray);

    /* create files */
    nc_file_define(file);

    double **f_x = dg_region_x(dg_phys_region(phys));
    double **f_y = dg_region_y(dg_phys_region(phys));

    /* set coordinate */
    nc_error( ncmpi_put_var_double_all(file->ncid, file->var_vec_p[0]->id, f_x[0]) );
    nc_error( ncmpi_put_var_double_all(file->ncid, file->var_vec_p[1]->id, f_y[0]) );

    extern Conv_Solver solver;
    solver.outfile = file;
    return;
}


void conv_putvar(dg_phys *phys, int timestep, double time){

    extern Conv_Solver solver;
    NC_File *file = solver.outfile;

    const int K = dg_grid_K(dg_phys_grid(phys));
    const int Np = dg_cell_Np(dg_phys_cell(phys));
    const int Nfield = dg_phys_Nfield(phys);

    dg_real *f_Q = dg_phys_f_Q(phys);

    MPI_Offset start_v[3], count_v[3];
    MPI_Offset start_t, count_t;

    NC_Var *ncvar = file->var_vec_p[2];
    /* put time */
    start_t = timestep; // start index
    count_t = 1;       // length
    nc_error( ncmpi_put_vara_double_all(file->ncid, ncvar->id, &start_t, &count_t, &time) );

    /* put variable c */
    start_v[0] = timestep;
    start_v[1] = 0;
    start_v[2] = 0;

    count_v[0] = 1;
    count_v[1] = K;
    count_v[2] = Np;

    ncvar = file->var_vec_p[3];
    float var[K*Np];
    /* put variable c */
    int k,n,sk=0,ind=0;
    for (k=0;k<K;k++){
        for (n=0;n<Np;n++){
            var[sk++] = (float) f_Q[ind]; // c field is the first
            ind += Nfield;
        }
    }
    nc_error( ncmpi_put_vara_float_all(file->ncid, ncvar->id, start_v, count_v, var) );
    return;
}


