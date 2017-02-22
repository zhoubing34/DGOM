//
// Created by li12242 on 17/2/20.
//

#include "pf_openbc.h"
#include "Utility/nc_library.h"
#include "Utility/utility.h"
#include "pf_phys.h"

/* global variable */
int name_flag=0;
char **filename = NULL;

/* local function */
static void time_linear_interp(double timeloc, int Nt, double *time, int *Ntstep, int *tstep, double *w);
static void pf_openbc_init(physField *phys, char *casename, int *ncfile);
static double* pf_openbc_interp(int ncid, double timeloc, time_interp_method method);

/**
 *
 * @param phys
 * @param casename
 * @param timeloc
 * @param method
 */
void pf_set_openbc(physField *phys, char *casename, double timeloc, time_interp_method method){
    parallMesh *mesh = phys->mesh;
    stdCell *cell = phys->cell;
    const int Nobc = mesh->Nobc;
    const int K = phys->grid->K;
    const int Nfield = phys->Nfield;
    const int Nv = cell->Nv;
    const int Np = cell->Np;
    int **EToV = phys->grid->EToV;
    real *f_ext = phys->f_ext;

    /* open boundary nc files */
    int *ncfile = vector_int_create(Nobc);
    pf_openbc_init(phys, casename, ncfile);

    int n,k,m,fld,sk=0;
    /* initialize of f_ext with the first obc */
    if(Nobc>0){
        double *bc = pf_openbc_interp(ncfile[0], timeloc, method);
        for(k=0;k<K;k++){
            for(fld=0;fld<Nfield;fld++){
                double vertbc[Nv];
                double nodebc[Np];
                for(m=0;m<Nv;m++){
                    int v1 = EToV[k][m];
                    vertbc[m] = bc[v1*Nfield + fld];
                }
                sc_vertProj(cell, vertbc, nodebc);

                sk = k*Np*Nfield +fld;
                for(m=0;m<Np;m++){
                    f_ext[sk + m*Nfield] = nodebc[m];
                }
            }
        }
        vector_double_free(bc);
        ncmpi_close(ncfile[0]);
    }
    /* loop over other obc */
    for(n=1;n<Nobc;n++){
        double *bc = pf_openbc_interp(ncfile[n], timeloc, method);
        // assign to f_ext
        for(k=0;k<K;k++){
            for(fld=0;fld<Nfield;fld++){
                double vertbc[Nv];
                double nodebc[Np];
                for(m=0;m<Nv;m++){
                    int v1 = EToV[k][m];
                    vertbc[m] = bc[v1*Nfield + fld];
                }
                sc_vertProj(cell, vertbc, nodebc);

                sk = k*Np*Nfield +fld;
                for(m=0;m<Np;m++){
                    f_ext[sk + m*Nfield] = nodebc[m];
                }
            }
        }
        vector_double_free(bc);
        ncmpi_close(ncfile[n]);
    }
    vector_int_free(ncfile);
    return;
}
/**
 * @brief create file names of open boundary
 * @param phys
 * @param casename
 */
static void pf_openbc_filename(physField *phys, char *casename){
    parallMesh *mesh = phys->mesh;
    const int Nobc = mesh->Nobc;
    int n;
    char indstr[20];
    filename = (char **)calloc(Nobc, sizeof(char*));
    for(n=0;n<Nobc;n++){
        filename[n] = (char *)calloc(MAX_NAME_LENGTH, sizeof(char));
        int ind = mesh->obcind[n];
        sprintf(indstr, "%d", ind);
        strcpy(filename[n], casename);
        strcat(filename[n], ".obc");
        strcat(filename[n], indstr);
        strcat(filename[n], ".nc");
    }
    name_flag = 1;
    return;
}
/**
 * @brief open boundary files
 * @param phys
 * @param casename
 * @param ncfile
 */
static void pf_openbc_init(physField *phys, char *casename, int *ncfile){
    parallMesh *mesh = phys->mesh;
    const int Nobc = mesh->Nobc;
    /* open boundary condition files */
    int n;
    for(n=0;n<Nobc;n++){
        if(!name_flag) pf_openbc_filename(phys, casename);
        int ret = ncmpi_open(MPI_COMM_WORLD, filename[n], NC_NOWRITE, MPI_INFO_NULL, ncfile+n);
        if(ret != NC_NOERR){
            fprintf(stderr, "%s(%s): %d\nError in open nc files: %s\n",
                    __FUNCTION__, __FILE__, __LINE__, ncmpi_strerror(ret));
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    return;
}

/**
 * @brief obtain the information of open boundary file
 * @param ncid
 * @param Nfield
 * @param Nt
 * @param time
 * @param Nv
 * @param vert
 */
void pf_openbc_info(int ncid, int *Nfield, int *Nt, double *time, int *Nv, int *vert){
    int varid, dimid;
    /* get time vector */
    MPI_Offset dimlen;
    ncmpi_inq_dimid(ncid, "Nt", &dimid); // get dimension id
    ncmpi_inq_dimlen(ncid, dimid, &dimlen); // get dimension length
    ncmpi_inq_varid(ncid, "time", &varid); // get variable id
    MPI_Offset start=0, count=dimlen;
    int nt = (int)dimlen;
    *Nt = nt;
    time = vector_double_create(nt);
    ncmpi_get_vara_double_all(ncid, varid, &start, &count, time);

    /* get vertex vector */
    ncmpi_inq_dimid(ncid, "Nv", &dimid); // get dimension id
    ncmpi_inq_dimlen(ncid, dimid, &dimlen); // get dimension length
    ncmpi_inq_varid(ncid, "vert", &varid); // get variable id
    const int nv = (int)dimlen;
    *Nv = nv;
    vert = vector_int_create(nv);
    start=0, count=dimlen;
    ncmpi_get_vara_int_all(ncid, varid, &start, &count, vert);

    /* get Nfield */
    ncmpi_inq_dimid(ncid, "Nfield", &dimid); // get dimension id
    ncmpi_inq_dimlen(ncid, dimid, &dimlen); // get dimension length
    *Nfield = (int)dimlen;
}

/**
 * @brief interp boundary values with given time
 * @param ncid
 * @param timeloc
 */
static double* pf_openbc_interp(int ncid, double timeloc, time_interp_method method){
    /* get basic info */
    int Nt, Nv, Nfield, *vert=NULL;
    double *time=NULL;
    pf_openbc_info(ncid, &Nfield, &Nt, time, &Nv, vert);

    /* determine time step and weights */
    int Nstep, *tstep=NULL;
    double *w=NULL;
    time_linear_interp(timeloc, Nt, time, &Nstep, tstep, w);

    /* get obc */
    int varid;
    ncmpi_inq_varid(ncid, "obc", &varid); // get variable id
    double *obc = vector_double_create(Nv*Nfield);
    double *tmp = vector_double_create(Nv*Nfield);
    MPI_Offset start_v[3] = {0, 0, 0};
    MPI_Offset count_v[3] = {1, Nv, Nfield};
    /* result from each time step */
    int t,n;
    for(t=0;t<Nstep;t++){
        start_v[0] = tstep[t];
        ncmpi_get_vara_double_all(ncid, varid, start_v, count_v, tmp);
        for(n=0;n<Nv*Nfield;n++)
            obc[n] += w[t]*tmp[n];
    }
    /* free memory */
    vector_double_free(time);
    vector_int_free(vert);
    vector_double_free(w);
    vector_int_free(tstep);
    vector_double_free(tmp);
    return obc;
}
/**
 * @brief adjacent time step and weights for interpolation
 * @param timeloc
 * @param Nt
 * @param time
 * @param Ntstep
 * @param tstep
 * @param w
 */
static void time_linear_interp(double timeloc, int Nt, double *time, int *Ntstep, int *tstep, double *w){
    const int nstep = 2;
    *Ntstep = nstep;
    w = vector_double_create(nstep);
    tstep = vector_int_create(nstep);
    /* determine time step and weights */
    if( (timeloc-time[0])<EPS ){
        tstep[0] = 0;
        tstep[1] = 0;
        w[0] = 0.5;
        w[1] = 0.5;
    }else if( (time[Nt-1]-timeloc)<EPS ){
        tstep[0] = Nt-1;
        tstep[1] = Nt-1;
        w[0] = 0.5;
        w[1] = 0.5;
    }else{
        int i;
        for(i=0;i<Nt;i++){
            if(timeloc < time[i]){
                tstep[1] = i;
                tstep[0] = i-1;
                double dt = time[i] - time[i-1];
                w[0] = (time[i] - timeloc)/dt;
                w[1] = (timeloc - time[i-1])/dt;
                break;
            }
        }
    }
    return;
}
