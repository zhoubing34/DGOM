#include "ConvectionDriver2d.h"
#include "NcOutput.h"

#ifndef DSET_NAME_LEN
#define DSET_NAME_LEN 1024
#endif
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
 * @warning
 * @attention
 * @note
 * @todo
 */
Ncfile* SetupOutput(MultiReg2d *mesh, char* casename){

    /* return file structure */
    Ncfile * outfile = (Ncfile *) calloc(1, sizeof(Ncfile));

    int ret, ncfile, ndims;
    StdRegions2d *shape = mesh->stdcell;

    /* dimension ids */
    int point_dimid, ele_dimid, time_dimid;
    int dimid2[2], dimid3[3];

    /* variable ids */
    int x_varid, y_varid, time_varid, tc_varid, var_varid;

    char filename[DSET_NAME_LEN];
    int procid = mesh->procid;
    int nprocs = mesh->nprocs;

#if defined DEBUG
    printf("Procs %d: Entering SetupOutput\n", procid);
#endif

    /* gen filename */
    ret = snprintf(filename, DSET_NAME_LEN, "%s.%d-%d.nc", casename, procid, nprocs);
    // check file name length
    if (ret >= DSET_NAME_LEN) {
        fprintf(stderr, "file name too long \n");
        exit(-1);
    }

    /* create output file */
    ret = ncmpi_create(MPI_COMM_SELF, filename,
                NC_CLOBBER|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncfile); NC_ERROR

    /* define dimensions */
    ret = ncmpi_def_dim(ncfile, "ele", mesh->K, &ele_dimid); NC_ERROR

    ret = ncmpi_def_dim(ncfile, "point", shape->Np, &point_dimid); NC_ERROR

    ret = ncmpi_def_dim(ncfile, "time", NC_UNLIMITED, &time_dimid); NC_ERROR

    /* define variables */
    dimid2[0] = ele_dimid;
    dimid2[1] = point_dimid;
    ndims = 2;
    ret = ncmpi_def_var(ncfile, "x", NC_DOUBLE, ndims, dimid2, &x_varid); NC_ERROR
    ret = ncmpi_def_var(ncfile, "y", NC_DOUBLE, ndims, dimid2, &y_varid); NC_ERROR

    ndims = 1;
    ret = ncmpi_def_var(ncfile, "time", NC_DOUBLE, ndims, &time_dimid, &time_varid); NC_ERROR

    /* set var */
    dimid3[0] = time_dimid;
    dimid3[1] = ele_dimid;
    dimid3[2] = point_dimid;
    ndims = 3;
    ret = ncmpi_def_var(ncfile, "var", NC_FLOAT, ndims, dimid3, &var_varid); NC_ERROR

    dimid2[0] = time_dimid;
    dimid2[1] = ele_dimid; // set var dimensions
    ndims = 2;
    ret = ncmpi_def_var(ncfile, "tcflag", NC_FLOAT, ndims, dimid2, &tc_varid); NC_ERROR

    /* end definations */
    ret = ncmpi_enddef(ncfile); NC_ERROR

    /* put node coordinate, x & y */
    ret = ncmpi_put_var_double_all(ncfile, x_varid, mesh->x[0]); NC_ERROR

    ret = ncmpi_put_var_double_all(ncfile, y_varid, mesh->y[0]); NC_ERROR

    /* set outfile file strucutre */
    outfile->ncfile = ncfile;
    outfile->varid = (int *)calloc(3, sizeof(int));

    outfile->varid[0] = time_varid;
    outfile->varid[1] = var_varid;
    outfile->varid[2] = tc_varid;

#if defined DEBUG
    printf("Procs %d: Leaving SetupOutput\n", procid);
#endif

    return outfile;
}


void PutVar(Ncfile * outfile, int outStep, double time, PhysDomain2d* phys){
    int ret;
    MPI_Offset start_f[2], count_f[2];
    MPI_Offset start_v[3], count_v[3];
    MPI_Offset start_t, count_t;
    MultiReg2d *mesh = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    /* put time */
    start_t = outStep; // start index
    count_t = 1;       // length
    ret = ncmpi_put_vara_double_all(outfile->ncfile, outfile->varid[0], &start_t, &count_t, &time); NC_ERROR

    /* put var */
    start_v[0] = outStep;
    start_v[1] = 0;
    start_v[2] = 0;
    count_v[0] = 1;
    count_v[1] = mesh->K;
    count_v[2] = shape->Np;
    ret = ncmpi_put_vara_float_all(outfile->ncfile, outfile->varid[1], start_v, count_v, phys->f_Q); NC_ERROR

//    /* put trouble cell flag */
//    start_f[0] = outStep;
//    start_f[1] = 0;
//    count_f[0] = 1;
//    count_f[1] = mesh->K;
//    ret = ncmpi_put_vara_float_all(outfile->ncfile, outfile->varid[2], start_v, count_v, mesh->tcflag); NC_ERROR

}


