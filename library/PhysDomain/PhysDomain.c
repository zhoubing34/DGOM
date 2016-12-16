/**
 * @file
 * Set physical variables
 *
 * @brief
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#include <MultiRegions/MultiRegBC/MultiRegBC2d.h>
#include "PhysDomain.h"

/* public variables */
void setSurfInfo(MultiRegBC2d *surf, int Nfields, real *surfinfo);
void permuteNodeIndexBuffer(MultiRegBC2d *, int Nfields, int *nodeIndexOut);
void permuteCellIndexBuffer(MultiRegBC2d *, int Nfields, int *cellIndexOut);

PhysDomain2d* PhysDomain2d_create(MultiReg2d *mesh, MultiRegBC2d *surf, int Nfields){

    PhysDomain2d *phys = (PhysDomain2d *) calloc(1, sizeof(PhysDomain2d));
    StdRegions2d *shape = mesh->stdcell;

    int K  = mesh->K;   /* number of elements */
    int Np = shape->Np; /* number of nodes in each element */

    int Nfp = shape->Nfp;
    int Nfaces = shape->Nfaces;

    /* number of variables */
    phys->Nfields = Nfields;

    /* mesh */
    phys->mesh = mesh;
    phys->surf = surf;

    /* volume info */
    phys->vgeo = mesh->vgeo;

    /* data array */
    phys->f_Q    = (real *) calloc(K*Np*Nfields, sizeof(real));
    phys->f_rhsQ = (real *) calloc(K*Np*Nfields, sizeof(real));
    phys->f_resQ = (real *) calloc(K*Np*Nfields, sizeof(real));
    phys->c_Q = (real *) calloc(K*Nfields, sizeof(real));
    /* external data */
    phys->f_ext = (real *) calloc(K*Np*Nfields, sizeof(real));

    /* node info send/recv buffer */
    int parNodeTotalOut = surf->parNodeTotalOut;
    phys->parNodeTotalOut = parNodeTotalOut*Nfields;
    phys->nodeIndexOut = IntVector_create(parNodeTotalOut*Nfields);
    size_t sz = (size_t) parNodeTotalOut*Nfields;
    phys->f_inQ  = (real *) calloc(sz, sizeof(real));
    phys->f_outQ = (real *) calloc(sz, sizeof(real));

    permuteNodeIndexBuffer(surf, Nfields, phys->nodeIndexOut);

    /* cell info send/recv buffer */
    int parCellTotalOut = surf->parCellTotalOut;
    phys->parCellTotalOut = parCellTotalOut*Nfields;
    phys->cellIndexOut = IntVector_create(parCellTotalOut*Nfields);
    sz = (size_t) parNodeTotalOut*Nfields;
    phys->c_inQ  = (real *) calloc(sz, sizeof(real));
    phys->c_outQ = (real *) calloc(sz, sizeof(real));

    permuteCellIndexBuffer(surf, Nfields, phys->cellIndexOut);

    /* surface info */
    int Nsurfinfo = 6;
    phys->Nsurfinfo = Nsurfinfo;
    sz = K*Nfp*Nfaces*Nsurfinfo*sizeof(real);
    phys->surfinfo = (real*) malloc(sz);
    setSurfInfo(surf, Nfields, phys->surfinfo);

    return phys;
}

void PhysDomain2d_free(PhysDomain2d *phys){

    free(phys->f_Q);
    free(phys->f_resQ);
    free(phys->f_rhsQ);
    free(phys->c_Q);
    free(phys->f_ext);

    free(phys->f_inQ);
    free(phys->f_outQ);

    free(phys->c_inQ);
    free(phys->c_outQ);

    IntVector_free(phys->nodeIndexOut);
    IntVector_free(phys->cellIndexOut);
    free(phys->surfinfo);
}

/**
 * @brief
 * Send and receive variables `f_Q` on boundaries in each processes.
 *
 * @details
 * The boundary values of `f_Q` on local process is arranged into `f_outQ` to send to other processes.
 * While `mpi_send_requests` get the request of each MPI_Isend function. The incoming boundary values
 * is stored into `f_inQ` and the request are stored in `mpi_recv_requests`.
 * `Nmess` gets the number of send and receive requests.
 *
 * @param[in]       phys    PhysDomain2d pointer
 * @param[inout]    mpi_send_requests  MPI_Request pointer
 * @param[inout]    mpi_send_requests  MPI_Request pointer
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * Nmessage | int      | number of messages stored in mpi_send_requests and mpi_send_requests
 *
 * Usages:
 *
 *     MPI_Request *mpi_out_requests = (MPI_Request*) calloc(mesh->nprocs, sizeof(MPI_Request));
 *     MPI_Request *mpi_in_requests  = (MPI_Request*) calloc(mesh->nprocs, sizeof(MPI_Request));
 *     int Nmess;
 *     fetchNodeBuffer2d(phys, mpi_out_requests, mpi_in_requests, &Nmess);
 *
 *     MPI_Status *instatus  = (MPI_Status*) calloc(nprocs, sizeof(MPI_Status));
 *     MPI_Waitall(Nmess, mpi_in_requests, instatus);
 *
 */
void fetchNodeBuffer2d(PhysDomain2d *phys,
                       MPI_Request *mpi_send_requests,
                       MPI_Request *mpi_recv_requests,
                       int *Nmessage) {
    int t;
    /* buffer outgoing node data */
    for(t=0;t<phys->parNodeTotalOut;++t)
        phys->f_outQ[t] = phys->f_Q[phys->nodeIndexOut[t]];

    MultiReg2d *mesh = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    /* do sends */
    int sk = 0, Nmess = 0;
    int p, Nout;
    for(p=0;p<mesh->nprocs;++p){
        if(p!=mesh->procid){
            Nout = mesh->Npar[p]*phys->Nfields*shape->Nfp; // # of variables send to process p
            if(Nout){
                /* symmetric communications (different ordering) */
                MPI_Isend(phys->f_outQ+sk, Nout, MPI_SIZE, p, 6666+p,
                          MPI_COMM_WORLD, mpi_send_requests +Nmess);
                MPI_Irecv(phys->f_inQ+sk,  Nout, MPI_SIZE, p, 6666+mesh->procid,
                          MPI_COMM_WORLD,  mpi_recv_requests +Nmess);
                sk+=Nout;
                ++Nmess;
            }
        }
    }
    *Nmessage = Nmess; /* number of messages */
}


void fetchCellBuffer(PhysDomain2d *phys,
                     MPI_Request *mpi_send_requests,
                     MPI_Request *mpi_recv_requests,
                     int *Nmessage){

    MultiReg2d *mesh = phys->mesh;

    /* buffer outgoing node data */
    int n;
    for(n=0;n<phys->parCellTotalOut;++n)
        phys->c_outQ[n] = phys->c_Q[phys->cellIndexOut[n]];

    /* do sends */
    int sk = 0, Nmess = 0, p;
    for(p=0;p<mesh->nprocs;++p){
        if(p!=mesh->procid){
            int Nout = mesh->Npar[p]; // # of variables send to process p
            if(Nout){
                /* symmetric communications (different ordering) */
                MPI_Isend(phys->c_outQ+sk, Nout, MPI_SIZE, p, 6666+p,
                          MPI_COMM_WORLD, mpi_send_requests +Nmess);
                MPI_Irecv(phys->c_inQ+sk,  Nout, MPI_SIZE, p, 6666+mesh->procid,
                          MPI_COMM_WORLD,  mpi_recv_requests +Nmess);
                sk+=Nout;
                ++Nmess;
            }
        }
    }
    *Nmessage = Nmess;
}

/**
 * @brief
 * Set the `nodeIndexOut` field for Phys2d structure.
 * @param[in] beginPos
 * @param[in] order order>0: year/month/date;order=0: date/month/year
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * car_id   | int      |
 * car_info | object   |
 *
 *
 */
void permuteNodeIndexBuffer(MultiRegBC2d *surf, int Nfields, int *nodeIndexOut){
    MultiReg2d *mesh = surf->mesh;
    int p2, n1, m, fld;
    int nprocs = mesh->nprocs;
    int procid = mesh->procid;
    StdRegions2d *shape = mesh->stdcell;
    int Nfp = shape->Nfp;

    int sp=0, sk=0;
    for(p2=0;p2<nprocs;++p2){
        if(p2!=procid) {
            /* for each received face */
            for (m=0;m<mesh->Npar[p2]; ++m) {
                for (n1=0; n1<Nfp; ++n1) {
                    for (fld = 0; fld < Nfields; ++fld) {
                        /* sk and node index determine the map relationship */
                        nodeIndexOut[sp++] = Nfields*(surf->nodeIndexOut[sk])+fld;
                    }
                    sk++;
                }
            }
        }
    }
}

void permuteCellIndexBuffer(MultiRegBC2d *surf, int Nfields, int *cellIndexOut){
    MultiReg2d *mesh = surf->mesh;
    int p2, m, fld;
    int nprocs = mesh->nprocs;
    int procid = mesh->procid;
    int sp=0, sk=0;
    for(p2=0;p2<nprocs;++p2){
        if(p2!=procid) {
            /* for each received face */
            for (m=0;m<mesh->Npar[p2];++m) {
                for (fld = 0; fld < Nfields;++fld) {
                    /* sk and node index determine the map relationship */
                    cellIndexOut[sp++] = Nfields*(surf->cellIndexOut[sk]) + fld;
                }
                sk++;
            }
        }
    }
}

void setSurfInfo(MultiRegBC2d *surf, int Nfields, real *surfinfo){

    int k,f,m,sk = 0;
    MultiReg2d *mesh = surf->mesh;
    StdRegions2d *shape = mesh->stdcell;

    int Np = shape->Np;
    int Nfaces = shape->Nfaces;
    int Nfp = shape->Nfp;

    double *drdx, *dsdx, *drdy, *dsdy, *J;
    int sz = Np*sizeof(double);

    drdx = (double*) malloc(sz);
    drdy = (double*) malloc(sz);
    dsdx = (double*) malloc(sz);
    dsdy = (double*) malloc(sz);
    J    = (double*) malloc(sz);

    /* local-local info */
    double *nxk = Vector_create(Nfaces);
    double *nyk = Vector_create(Nfaces);
    double *sJk = Vector_create(Nfaces);

    for(k=0;k<mesh->K;++k){
        GeoFactor2d(shape->Np, mesh->x[k], mesh->y[k],
                    shape->Dr, shape->Ds,
                    drdx, dsdx, drdy, dsdy, J);
        Normals2d(Nfaces, mesh->GX[k], mesh->GY[k], nxk, nyk, sJk);

        for(f=0;f<Nfaces;++f){
            int bcType = surf->EToBS[k][f];
            for(m=0;m<Nfp;++m){

                int id  = m + f*Nfp + Nfp*Nfaces*k;
                int idM = surf->vmapM[id];
                int idP = surf->vmapP[id];
                int  nM = idM%Np;
                int  nP = idP%Np;
                int  kM = (idM-nM)/Np;
                int  kP = (idP-nP)/Np;

                idM = Nfields*(nM+Np*kM);
                idP = Nfields*(nP+Np*kP);

                /* stub resolve some other way */
                if(surf->vmapP[id]<0){
                    idP = surf->vmapP[id]; /* -ve numbers */
                    idP = Nfields*(-1-idP);
                }

                surfinfo[sk++] = idM;
                surfinfo[sk++] = idP;
                surfinfo[sk++] = (real)(sJk[f]/(J[shape->Fmask[f][m]]));
                surfinfo[sk++] = bcType;
                surfinfo[sk++] = (real)nxk[f];
                surfinfo[sk++] = (real)nyk[f];
            }
        }
    }
    /* deallocate mem */
    free(drdx);
    free(drdy);
    free(dsdx);
    free(dsdy);
    free(J);

    Vector_free(nxk);
    Vector_free(nyk);
    Vector_free(sJk);

}
