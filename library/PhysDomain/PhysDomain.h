#ifndef DGOM_PHYSDOMAIN_H
#define DGOM_PHYSDOMAIN_H

#include "MultiRegions/MultiRegions.h"

typedef struct PhysDomain2d{
    /** number of variable fields */
    int Nfields;
    /** mesh */
    MultiReg2d* mesh;
    /** surf infomation */
    real *surfinfo;
    /** volume info */
    real *vgeo;

    /** number of variables to send/recv */
    int parNtotalout;
    /** map from send/recv array to variables */
    int *parmapOUT;
    /** variable array */
    real *f_Q, *f_rhsQ, *f_resQ;

    /** send/recv data */
    real *f_inQ, *f_outQ;
}PhysDomain2d;

PhysDomain2d* GetPhysDomain2d(MultiReg2d *mesh, int Nfields);
void FreePhysDomain2d(PhysDomain2d *phys);

void FetchParmapNode2d(PhysDomain2d *phys,
                       MPI_Request *mpi_send_requests,
                       MPI_Request *mpi_recv_requests,
                       int *Nmessage);

void FetchParmapEle2d(PhysDomain2d *phys, float *f_E,
                      float *f_inE, float *f_outE,
                      MPI_Request *mpi_send_requests,
                      MPI_Request *mpi_recv_requests,
                      int *Nmessage);

#endif //DGOM_PHYSDOMAIN_H
