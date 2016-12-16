#ifndef DGOM_PHYSDOMAIN_H
#define DGOM_PHYSDOMAIN_H

#include "MultiRegions/MultiRegions.h"
#include "MultiRegions/MultiRegBC/MultiRegBC2d.h"

typedef struct{
    /** number of variable fields */
    int Nfields;
    /** mesh */
    MultiReg2d *mesh;
    MultiRegBC2d *surf;

    /** number of surf infomation */
    int Nsurfinfo;
    /** surf infomation */
    real *surfinfo;
    /** number of volume information */
    int Nvgeo;
    /** volume info */
    real *vgeo;

    /** number of node variables to send/recv */
    int parNodeTotalOut;
    /** map from send/recv array to variables */
    int *nodeIndexOut;
    /** send/recv data */
    real *f_inQ, *f_outQ;

    /** number of cell variables to send/recv */
    int parCellTotalOut;
    /** array of cell index to send/recv  */
    int *cellIndexOut;
    /** send/recv data */
    real *c_inQ, *c_outQ;

    /** volume info */
    real *c_Q;

    /** external data on vertex */
    real *f_ext;

    /** variable array */
    real *f_Q, *f_rhsQ, *f_resQ;
}PhysDomain2d;

PhysDomain2d* PhysDomain2d_create(MultiReg2d *mesh, MultiRegBC2d *surf, int Nfields);
void PhysDomain2d_free(PhysDomain2d *phys);

void fetchNodeBuffer2d(PhysDomain2d *phys,
                       MPI_Request *mpi_send_requests,
                       MPI_Request *mpi_recv_requests,
                       int *Nmessage);

void fetchCellBuffer(PhysDomain2d *phys,
                     MPI_Request *mpi_send_requests,
                     MPI_Request *mpi_recv_requests,
                     int *Nmessage);

#endif //DGOM_PHYSDOMAIN_H
