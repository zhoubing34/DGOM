//
// Created by li12242 on 17/2/20.
//

#include "dg_phys_obc.h"

#define DEBUG 1

typedef struct Interp_Parameter{
    int Nstep; ///< number of time step;
    int *timeStep; ///< index of time step;
    double *weight; ///< interp parameters;
}Interp_Parameter;

static void dg_phys_obc_add(dg_phys_obc *phys_obc, char *filename);
static void dg_phys_obc_linear_interp(dg_phys_obc *phys_obc, Interp_Parameter *par, double elapseTime);
static void dg_phys_obc_obtain(dg_phys_obc *phys_obc, Interp_Parameter *par);
static void dg_phys_obc_update(dg_phys_obc *phys_obc, double elapseTime);

dg_phys_obc* dg_phys_obc_create(dg_phys_info *phys_info){
    dg_phys_obc *phys_obc = calloc(1, sizeof(dg_phys_obc));

    phys_obc->info = phys_info;
    phys_obc->f_extQ = phys_info->f_Q; // (default)
    phys_obc->interp_type = INTERP_LINEAR; // (default)
    phys_obc->file = NULL; // (default)

    phys_obc->add_obc = dg_phys_obc_add;
    phys_obc->update_obc = dg_phys_obc_update;
    return phys_obc;
}

void dg_phys_obc_free(dg_phys_obc *phys_obc){

    if(phys_obc->file != NULL) {
        nc_file_free(phys_obc->file);
        free(phys_obc->f_extQ);
    }
    free(phys_obc);
    return;
}

/**
 *
 * @param phys
 * @param casename
 * @param timeloc
 * @param method
 */
static void dg_phys_obc_add(dg_phys_obc *phys_obc, char *filename){

    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    NC_File *file = nc_file_read_from_file(filename, procid, nprocs);

    /* check obc file variables */
    NC_Var *varVert = file->var_vec_p[0];
    NC_Var *varTime = file->var_vec_p[1];
    NC_Var *varOBC = file->var_vec_p[2];
    if( memcmp(varVert->name, "vert", sizeof("vert")) )
        fprintf(stderr, "%s (%d): The first variable in open boundary file "
                "(%s) shoule be 'vert'\n", __FUNCTION__, __LINE__, filename);
    if( memcmp(varTime->name, "time", sizeof("time")) )
        fprintf(stderr, "%s (%d): The first variable in open boundary file "
                "(%s) shoule be 'time'\n", __FUNCTION__, __LINE__, filename);
    if( memcmp(varOBC->name, "obc", sizeof("obc")) )
        fprintf(stderr, "%s (%d): The first variable in open boundary file "
                "(%s) shoule be 'obc'\n", __FUNCTION__, __LINE__, filename);
    /* check input Nfield */
    int Nfield = varOBC->dim_vec_p[2]->len;
    if(Nfield != phys_obc->info->Nfield)
        fprintf(stderr, "%s (%d): The Nfield (%d) in obc file is incorrect\n",
                __FUNCTION__, __LINE__, Nfield);

    /* variable information */
    int Nt = varOBC->dim_vec_p[0]->len;
    int Nvert = varOBC->dim_vec_p[1]->len;

    int ncid = file->ncid;
    MPI_Offset start=0, count=Nt;
    double *time = vector_double_create(Nt);
    int *vert = vector_int_create(Nvert);
    /* read time */
    ncmpi_get_vara_double_all(ncid, varTime->id, &start, &count, time);
    /* read vertex list */
    start=0, count=Nvert;
    ncmpi_get_vara_int_all(ncid, varVert->id, &start, &count, vert);

    /* assignment */
    phys_obc->Ntime = Nt;
    phys_obc->Nvert = Nvert;
    phys_obc->time = time;
    phys_obc->vert = vert;
    phys_obc->file = file;

    /* allocation */
    const int K = dg_grid_K(dg_phys_info_grid(phys_obc->info));
    const int Np = dg_cell_Np(dg_phys_info_cell(phys_obc->info));
    phys_obc->f_extQ = (dg_real *) calloc((size_t) K*Np*Nfield, sizeof(dg_real));
    return;
}

/**
 * @brief
 * @param phys_obc
 * @param elapseTime
 */
static void dg_phys_obc_update(dg_phys_obc *phys_obc, double elapseTime){
    /* check file */
    if(phys_obc->file == NULL){
        return;
    }

    Interp_Parameter *par = calloc(1, sizeof(Interp_Parameter));
    switch (phys_obc->interp_type){
        case INTERP_LINEAR:
            dg_phys_obc_linear_interp(phys_obc, par, elapseTime);
            break;
        default:
            fprintf(stderr, "%s (%d): Unknown interpolation type\n", __FUNCTION__, __LINE__);
            break;
    }
    dg_phys_obc_obtain(phys_obc, par);
    free(par);
    return;
}

static void dg_phys_obc_obtain(dg_phys_obc *phys_obc, Interp_Parameter *par){

    const int ncid = phys_obc->file->ncid;
    const int varid = phys_obc->file->var_vec_p[2]->id;
    const int Nfield = phys_obc->info->Nfield;
    const int Nvobc = phys_obc->Nvert; ///< number of vertex in obc file;
    const int Nvert = dg_grid_Nv(dg_phys_info_grid(phys_obc->info)); ///< number of vertex on geometry grid;
    const int Nstep = par->Nstep;
    double *weight = par->weight;
    int *tstep = par->timeStep;
    int *vert = phys_obc->vert;
    /* get open boundary vertex result */
    double *obc = vector_double_create(Nvert*Nfield);
    double *tmp = vector_double_create(Nvobc*Nfield);
    MPI_Offset start_v[3] = {0, 0, 0};
    MPI_Offset count_v[3] = {1, Nvobc, Nfield};
    /* result from each time step */
    int t,n,k,m,fld;
    for(t=0;t<Nstep;t++){
        double wtmp = weight[t];
        start_v[0] = tstep[t];
        ncmpi_get_vara_double_all(ncid, varid, start_v, count_v, tmp);
        for(n=0;n<Nvobc;n++){
            int sk = vert[n]*Nfield;
            int sp = n*Nfield;
            for(fld=0;fld<Nfield;fld++){
                obc[sk+fld] += wtmp*tmp[sp+fld];
            }
        }
    }
    /* interp the open boundary vertex onto nodes */
    dg_cell *cell = dg_phys_info_cell(phys_obc->info);
    const int K = dg_grid_K(dg_phys_info_grid(phys_obc->info));
    const int Nv = dg_cell_Nv(cell);
    const int Np = dg_cell_Np(cell);

    double *f_ext = phys_obc->f_extQ;
    int **EToV = dg_grid_EToV(dg_phys_info_grid(phys_obc->info));
    for(k=0;k<K;k++){
        double vertbc[Nv*Nfield];
        double nodebc[Np*Nfield];
        for(m=0;m<Nv;m++){
            int v1 = EToV[k][m]; // vertex index
            for(fld=0;fld<Nfield;fld++){
                vertbc[m*Nfield+fld] = obc[v1*Nfield+fld];
            }
        }
        cell->proj_vert2node(cell, Nfield, vertbc, nodebc);
        // assignment to nodes value
        int sk = k*Np*Nfield;
        for(m=0;m<Np*Nfield;m++){
            f_ext[sk+m] = (dg_real) nodebc[m];
        }
    }
    return;
}

/**
 * @brief adjacent time step and weights for interpolation
 * @param elapseTime
 * @param Nt
 * @param time
 * @param Ntstep
 * @param tstep
 * @param w
 */
static void dg_phys_obc_linear_interp(dg_phys_obc *phys_obc, Interp_Parameter *par, double elapseTime){

    const int Nstep = 2;
    const int Nt = phys_obc->Ntime;

    double *w = vector_double_create(Nstep);
    int *tstep = vector_int_create(Nstep);
    int n;
    double *time = phys_obc->time;
    /* determine time step and weights */
    if( (elapseTime-time[0])<EPS ){
        tstep[0] = 0;
        tstep[1] = 0;
        w[0] = 0.5;
        w[1] = 0.5;
    }else if( (time[Nt-1]-elapseTime)<EPS ){
        tstep[0] = Nt-1;
        tstep[1] = Nt-1;
        w[0] = 0.5;
        w[1] = 0.5;
    }else{
        for(n=0;n<Nt;n++){
            if(elapseTime < time[n]){
                tstep[1] = n;
                tstep[0] = n-1;
                double dt = time[n] - time[n-1];
                w[0] = (time[n] - elapseTime)/dt;
                w[1] = (elapseTime - time[n-1])/dt;
                break;
            }
        }
    }
    // assignment
    par->Nstep = Nstep;
    par->timeStep = tstep;
    par->weight = w;
    return;
}