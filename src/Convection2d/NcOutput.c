#include "Convection2d/Convection2d.h"

#define DSET_NAME_LEN 1024

/*
 * create output files in NetCDF format
 * one file per process
 * */
Ncfile* SetupOutput(Mesh *mesh, char* filename){
    Ncfile * outfile = (Ncfile *) calloc(1, sizeof(Ncfile));

    int ret, ncfile, ndims;
    int node_dim, time_dim, var_dim[2];
    int xid, yid, timeid, varid;

    ret = snprintf(filename, DSET_NAME_LEN, "%s.%d-%d.nc", filename, mesh->procid, mesh->nprocs);
    // check file name length
    if (ret >= DSET_NAME_LEN) {
        fprintf(stderr, "file name too long \n");
        exit(-1);
    }

    // create output file
    ret = ncmpi_create(MPI_COMM_SELF, filename,
                NC_CLOBBER|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncfile);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    // define dimensions
    ret = ncmpi_def_dim(ncfile, "node", p_Np*mesh->K, &node_dim);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    ret = ncmpi_def_dim(ncfile, "time", NC_UNLIMITED, &time_dim);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    // define variables
    ndims = 1;
    ret = ncmpi_def_var(ncfile, "x", NC_DOUBLE, ndims, &node_dim, &xid);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    ret = ncmpi_def_var(ncfile, "y", NC_DOUBLE, ndims, &node_dim, &yid);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    ret = ncmpi_def_var(ncfile, "time", NC_DOUBLE, ndims, &time_dim, &timeid);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    var_dim[0] = time_dim; var_dim[1] = node_dim; // set var dimensions
    ndims = 2;
    ret = ncmpi_def_var(ncfile, "var", NC_FLOAT, ndims, var_dim, &varid);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    // end definations
    ret = ncmpi_enddef(ncfile); if (ret != NC_NOERR) handle_error(ret, __LINE__);

    // put node coordinate, x & y
    ret = ncmpi_put_var_double_all(ncfile, xid, mesh->x[0]);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    ret = ncmpi_put_var_double_all(ncfile, yid, mesh->y[0]);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    // set to outfile field
    outfile->ncfile = ncfile;
    outfile->varid = (int *)calloc(2, sizeof(int));

    outfile->varid[0] = timeid; outfile->varid[1] = varid;

    return outfile;
}

void PutVar(Ncfile * outfile, int outStep, double time, Mesh* mesh){
    int ret;
    MPI_Offset start_v[2], count_v[2];
    MPI_Offset start_t, count_t;

    // put time
    start_t = outStep; // start index
    count_t = 1;       // length
    ret = ncmpi_put_vara_double_all(outfile->ncfile, outfile->varid[0], &start_t, &count_t, &time);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    // put var
    start_v[0] = outStep;
    start_v[1] = 0;
    count_v[0] = 1;
    count_v[1] = p_Np*mesh->K;
    ret = ncmpi_put_vara_float_all(outfile->ncfile, outfile->varid[1], start_v, count_v, mesh->f_Q);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

}


