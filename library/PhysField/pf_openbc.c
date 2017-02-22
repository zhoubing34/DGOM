//
// Created by li12242 on 17/2/20.
//

#include "pf_openbc.h"
#include "Utility/nc_library.h"

void pf_get_openbc(physField *phys, double timeloc, time_interp_method method){
    parallMesh *mesh = phys->mesh;
    /* open boundary condition files */

}

static void linear_interp(int ncid, double tloc){

    /* get time variable */
    int varid, dimid;
    MPI_Offset dimlen;
    ncmpi_inq_dimid(ncid, "Nt", &dimid); // get dimension id
    ncmpi_inq_dimlen(ncid, dimid, &dimlen); // get dimension length
    ncmpi_inq_varid(ncid, "time", &varid); // get variable id
    MPI_Offset start=0, count=dimlen;
    const int Nt = (int)dimlen;
    double *time = vector_double_create(Nt);
    ncmpi_get_vara_double_all(ncid, varid, &start, &count, time);

    /* get vertex list variable */
    ncmpi_inq_dimid(ncid, "Nv", &dimid); // get dimension id
    ncmpi_inq_dimlen(ncid, dimid, &dimlen); // get dimension length
    ncmpi_inq_varid(ncid, "vert", &varid); // get variable id
    const int Nv = (int)dimlen;
    int *vert = vector_int_create(Nv);
    start=0, count=dimlen;
    ncmpi_get_vara_int_all(ncid, varid, &start, &count, vert);

    ncmpi_inq_dimid(ncid, "Nfield", &dimid); // get dimension id
    ncmpi_inq_dimlen(ncid, dimid, &dimlen); // get dimension length
    const int Nfield = (int)dimlen;

    /* determine time step */
    int i, pre_tstep=0, next_tstep=0;
    double w_pre, w_next;
    if( (tloc-time[0])<EPS ){
        next_tstep = 0;
        pre_tstep = 0;
        w_pre = 0.5;
        w_next = 0.5;
    }else if( (time[Nt-1]-tloc)<EPS ){
        next_tstep = Nt-1;
        pre_tstep = Nt-1;
        w_pre = 0.5;
        w_next = 0.5;
    }else{
        for(i=0;i<Nt;i++){
            if(tloc < time[i]){
                next_tstep = i;
                pre_tstep = i-1;
                double dt = time[i] - time[i-1];
                w_pre = (time[i] - tloc)/dt;
                w_next = (tloc - time[i-1])/dt;
                break;
            }
        }
    }

    /* get obc */
    ncmpi_inq_varid(ncid, "obc", &varid); // get variable id
    double *obc = vector_double_create(Nv*Nfield);
    MPI_Offset start_v[3], count_v[3];
    start_v[0] = pre_tstep;
    start_v[1] = 0;
    start_v[2] = 0;
    count_v[0] = 1;
    count_v[1] = Nv;
    count_v[2] = Nfield;
    ncmpi_get_vara_double_all(ncid, varid, start_v, count_v, obc);

    vector_double_free(time);
    vector_int_free(vert);
    vector_double_free(obc);
    return;
}
