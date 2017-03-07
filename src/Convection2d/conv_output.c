#include <MultiRegions/Mesh/mr_mesh.h>
#include <PhysField/pf_phys.h>
#include "conv_driver2d.h"

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
void conv_setoutput(physField *phys){

    const int K = phys->grid->K;
    const int Np = phys->cell->Np;
    const int procid = phys->mesh->procid;
    const int nprocs = phys->mesh->nprocs;

    /* define dimensions */
    nc_dim *ne = nc_dim_create("ne", K);
    nc_dim *np = nc_dim_create("np", Np);
    nc_dim *t  = nc_dim_create("t", 0);

    const int ndim = 3;
    const int nvar = 4;
    nc_dim *dimarray[ndim];
    nc_var *vararray[nvar];

    /* define variables: x,y,time and h */
    dimarray[0] = ne;
    dimarray[1] = np; /* the inner loop dimension comes the last */
    nc_var *x = nc_var_create("x", 2, dimarray, NC_DOUBLE);
    nc_var *y = nc_var_create("y", 2, dimarray, NC_DOUBLE);

    dimarray[0] = t;
    nc_var *time = nc_var_create("time", 1, dimarray, NC_DOUBLE);

    dimarray[0] = t;
    dimarray[1] = ne;
    dimarray[2] = np; /* inner loop */
    nc_var *h = nc_var_create("var", 3, dimarray, NC_FLOAT);

    /* define files */
    dimarray[0] = np;
    dimarray[1] = ne;
    dimarray[2] = t;

    vararray[0] = x;
    vararray[1] = y;
    vararray[2] = time;
    vararray[3] = h;
    nc_file *file = nc_file_create("conv2d.", procid, nprocs, ndim, dimarray, nvar, vararray);

    /* create files */
    nc_file_init(file);

    /* set coordinate */
    nc_error( ncmpi_put_var_double_all(file->id, file->var_vec_p[0]->id, phys->region->x[0]) );
    nc_error( ncmpi_put_var_double_all(file->id, file->var_vec_p[1]->id, phys->region->y[0]) );

    extern conv_solver2d solver;
    solver.outfile = file;

    return;
}


void conv_putvar(physField *phys, int timestep, double time){

    extern conv_solver2d solver;
    nc_file *file = solver.outfile;

    const int K = phys->grid->K;
    const int Np = phys->cell->Np;
    const int Nfield = phys->Nfield;

    MPI_Offset start_v[3], count_v[3];
    MPI_Offset start_t, count_t;

    nc_var *ncvar = file->var_vec_p[2];
    /* put time */
    start_t = timestep; // start index
    count_t = 1;       // length
    nc_error( ncmpi_put_vara_double_all(file->id, ncvar->id, &start_t, &count_t, &time) );

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
            var[sk++] = (float) phys->f_Q[ind]; // c field is the first
            ind += Nfield;
        }
    }
    nc_error( ncmpi_put_vara_float_all(file->id, ncvar->id, start_v, count_v, var) );

    return;
}


