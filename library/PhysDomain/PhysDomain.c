/**
 * @file
 * Set physical variables
 *
 * @brief
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "PhysDomain.h"

/* public variables */
void SetSurfInfo2d(MultiReg2d *mesh, int Nfields, real *surfinfo);
void SetParmapOut2d(MultiReg2d *mesh, int Nfields, int *parmapOUT);

PhysDomain2d* GenPhysDomain2d(MultiReg2d *mesh, int Nfields){

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

    /* volume info */
    phys->vgeo = mesh->vgeo;

    /* data array */
    phys->f_Q    = (real *) calloc(K*Np*Nfields, sizeof(real));
    phys->f_rhsQ = (real *) calloc(K*Np*Nfields, sizeof(real));
    phys->f_resQ = (real *) calloc(K*Np*Nfields, sizeof(real));

    /* MPI send/recv buffer */
    phys->f_inQ  = (real *) calloc(mesh->parNtotalout*Nfields, sizeof(real));
    phys->f_outQ = (real *) calloc(mesh->parNtotalout*Nfields, sizeof(real));

    phys->parNtotalout = mesh->parNtotalout*Nfields;
    phys->parmapOUT = IntVector_create(phys->parNtotalout);
    SetParmapOut2d(mesh, Nfields, phys->parmapOUT);

    /* surface info */
    int sz = K*Nfp*Nfaces*6*sizeof(real);
    phys->surfinfo = (real*) malloc(sz);
    SetSurfInfo2d(mesh, Nfields, phys->surfinfo);

    return phys;
}

void FreePhysDomain2d(PhysDomain2d *phys){
    free(phys->f_Q);
    free(phys->f_resQ);
    free(phys->f_rhsQ);
    free(phys->f_inQ);
    free(phys->f_outQ);

    IntVector_free(phys->parmapOUT);
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
 *     FetchParmapNode2d(phys, mpi_out_requests, mpi_in_requests, &Nmess);
 *
 *     MPI_Status *instatus  = (MPI_Status*) calloc(nprocs, sizeof(MPI_Status));
 *     MPI_Waitall(Nmess, mpi_in_requests, instatus);
 *
 */
void FetchParmapNode2d(PhysDomain2d *phys,
                       MPI_Request *mpi_send_requests,
                       MPI_Request *mpi_recv_requests,
                       int *Nmessage) {
    int t;
    /* buffer outgoing node data */
    for(t=0;t<phys->parNtotalout;++t)
        phys->f_outQ[t] = phys->f_Q[phys->parmapOUT[t]];

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


void FetchParmapEle2d(PhysDomain2d *phys,
                      real *f_E, real *f_inE, real *f_outE,
                      MPI_Request *mpi_send_requests,
                      MPI_Request *mpi_recv_requests,
                      int *Nmessage){

    MultiReg2d *mesh = phys->mesh;

    /* buffer outgoing node data */
    int n;
    for(n=0;n<mesh->parEtotalout;++n)
        f_outE[n] = f_E[mesh->elemapOut[n]];

    /* do sends */
    int sk = 0, Nmess = 0, p;
    for(p=0;p<mesh->nprocs;++p){
        if(p!=mesh->procid){
            int Nout = mesh->Npar[p]; // # of variables send to process p
            if(Nout){
                /* symmetric communications (different ordering) */
                MPI_Isend(f_outE+sk, Nout, MPI_SIZE, p, 6666+p,
                          MPI_COMM_WORLD, mpi_send_requests +Nmess);
                MPI_Irecv(f_inE+sk,  Nout, MPI_SIZE, p, 6666+mesh->procid,
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
 * Set the `parmapOUT` field for PhysDomain2d structure.
 *
 * @details
 *
 *
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
void SetParmapOut2d(MultiReg2d *mesh, int Nfields, int *parmapOUT){
    int p2, n1, m, fld;
    int nprocs = mesh->nprocs;
    int procid = mesh->procid;
    StdRegions2d *shape = mesh->stdcell;
    int Nfp = shape->Nfp;

    int sp=0, sk=0;
    for(p2=0;p2<nprocs;++p2){
        if(p2!=procid) {
            /* for each received face */
            for (m = 0; m < mesh->Npar[p2]; ++m) {
                for (n1 = 0; n1 < Nfp; ++n1) {
                    for (fld = 0; fld < Nfields; ++fld) {
                        /* sk and node index determine the map relationship */
                        parmapOUT[sp++] = Nfields * (mesh->parmapOUT[sk]) + fld;
                    }
                    sk++;
                }
            }
        }
    }
}


void SetSurfInfo2d(MultiReg2d *mesh, int Nfields, real *surfinfo){

    int k,f,m,sk = 0;
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
            for(m=0;m<Nfp;++m){

                int id  = m + f*Nfp + Nfp*Nfaces*k;
                int idM = mesh->vmapM[id];
                int idP = mesh->vmapP[id];
                int  nM = idM%Np;
                int  nP = idP%Np;
                int  kM = (idM-nM)/Np;
                int  kP = (idP-nP)/Np;

                idM = Nfields*(nM+Np*kM);
                idP = Nfields*(nP+Np*kP);

                /* stub resolve some other way */
                if(mesh->vmapP[id]<0){
                    idP = mesh->vmapP[id]; /* -ve numbers */
                }

                surfinfo[sk++] = idM;
                surfinfo[sk++] = idP;
                surfinfo[sk++] = (real)(sJk[f]/(J[shape->Fmask[f][m]]));
                surfinfo[sk++] = (real)sJk[f];
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
