#include "Convection2d/Convection2d.h"
#include <mpi.h>

void ConvectionRHS2d(Mesh *mesh, float frka, float frkb, float fdt){
/* registers and temporary */
    register unsigned int k, n;

    /* mesh parameters */
    const int K = mesh->K;

    float *vgeo     = mesh->vgeo;
    float *surfinfo = mesh->surfinfo;
    float *f_Dr    = mesh->f_Dr;
    float *f_Ds    = mesh->f_Ds;
    float *f_LIFT  = mesh->f_LIFT;

    float *f_Q    = mesh->f_Q;
    float *f_rhsQ = mesh->f_rhsQ;
    float *f_resQ = mesh->f_resQ;
    float *f_s    = mesh->f_s;

    float *f_inQ  = mesh->f_inQ;
    float *f_outQ  = mesh->f_outQ;

    int p;

    /* mpi request buffer */
    MPI_Request *mpi_out_requests = (MPI_Request*) calloc(mesh->nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_in_requests  = (MPI_Request*) calloc(mesh->nprocs, sizeof(MPI_Request));

    /* buffer outgoing node data */
    for(n=0;n<mesh->parNtotalout;++n)
        mesh->f_outQ[n] = f_Q[mesh->parmapOUT[n]];

    /* do sends */
    int sk = 0, Nmess = 0, uk = 0;
    for(p=0;p<mesh->nprocs;++p){
        if(p!=mesh->procid){
            int Nout = mesh->Npar[p]*p_Nfields*p_Nfp; // # of variables send to process p
            if(Nout){
                /* symmetric communications (different ordering) */
                MPI_Isend(f_outQ+sk, Nout, MPI_FLOAT, p, 6666+p,            MPI_COMM_WORLD, mpi_out_requests +Nmess);
                MPI_Irecv(f_inQ+sk,  Nout, MPI_FLOAT, p, 6666+mesh->procid, MPI_COMM_WORLD,  mpi_in_requests +Nmess);
                sk+=Nout;
                ++Nmess;
            }// if
        }// if
    }// for

    // volume integral
    for(k=0;k<K;++k){

        /* NOTE: once k is known, all other indexing variables etc are derived */
        register unsigned int n, m;

        /* NOTE: should be local memory */
        float  Qk[p_Np*p_Nfields];
        float  Uf[p_Np*2];

        /* NOTE: index into geometric factors */
        int geoid=k*4;

        const float drdx = vgeo[geoid++], drdy = vgeo[geoid++];
        const float dsdx = vgeo[geoid++], dsdy = vgeo[geoid++];

        /* NOTE: buffer element k into local storage */
        float *qpt = f_Q + p_Nfields*p_Np*k;
        float *upt = f_s + 2*p_Np*k;

        uk = 0;
        for(m=0;m<p_Np;++m){
            Uf[uk++] = upt[uk];
            Uf[uk++] = upt[uk];
        }

        for(m=0;m<p_Nfields*p_Np;++m){
            Qk[m] = qpt[m];
        }

        for(n=0;n<p_Np;++n){

            const float *ptDr = f_Dr+n*p_Np;
            const float *ptDs = f_Ds+n*p_Np;

            float rhs = 0;

            sk = 0; uk = 0;
            for(m=0;m<p_Np;++m){
                const float dr = ptDr[m];
                const float ds = ptDs[m];
                const float dx = drdx*dr+dsdx*ds;
                const float dy = drdy*dr+dsdy*ds;
                const float u = Uf[uk++];
                const float v = Uf[uk++];

                const float C = Qk[sk++];

                rhs += - dx*u*C;
                rhs += - dy*v*C;
            }

            int id = p_Nfields*(k*p_Np + n);
            f_rhsQ[id] = rhs;
        }// for n
    }// for k

    /* DO RECV */
    MPI_Status *instatus  = (MPI_Status*) calloc(mesh->nprocs, sizeof(MPI_Status));
    MPI_Waitall(Nmess, mpi_in_requests, instatus);
    free(instatus);

    // surface integral
    for(k=0;k<K;++k){

        /* NOTE: once k is known, all other indexing variables etc are derived */
        register unsigned int n, m;

        /* NOTE: should be local memory */
        float fluxQ[p_Nfaces*p_Nfp*p_Nfields];

        /* NOTE: index into geometric factors */
        int surfid=k*6*p_Nfp*p_Nfaces;

        /* Lax-Friedrichs flux */
        int sk = 0;
        for(m=0;m<p_Nfp*p_Nfaces;++m){

            int   idM       = surfinfo[surfid++];
            int   idP       = surfinfo[surfid++];
            const float FSc = surfinfo[surfid++];
            const float BSc = surfinfo[surfid++];
            const float NXf = surfinfo[surfid++];
            const float NYf = surfinfo[surfid++];

            float dC, dF, dG;
            if(idP<0){
                idP = p_Nfields*(-1-idP);
                dC = (f_Q[idP] - f_inQ[idM]);
            }
            else{
                dC = ( f_Q[idP] - f_Q[idM] );
            }

            const float uM = f_s[idM*2];
            const float vM = f_s[idM*2+1];

            dF = uM*dC;
            dG = vM*dC;

            const float un = fabsf(NXf*uM + NYf*vM);
            fluxQ[sk++] = FSc*( - NXf*dF - NYf*dG + un*dC);

        }

        /* Lift the flux term */
        for(n=0;n<p_Np;++n){

            const float *ptLIFT = f_LIFT+n*p_Nfp*p_Nfaces;

            float rhs = 0;

            sk = 0;
            for(m=0;m<p_Nfp*p_Nfaces;++m){
                const float L = ptLIFT[m];
                rhs += L*fluxQ[sk++];
            }

            int id = p_Nfields*(n+k*p_Np);
            f_rhsQ[id] += rhs;
        }
    }

    for(n=0;n<K*p_Np*p_Nfields;++n){
        f_resQ[n] = frka*f_resQ[n]+fdt*f_rhsQ[n];
        f_Q[n]   += frkb*f_resQ[n];
    }

    /* make sure all messages went out */
    MPI_Status *outstatus  = (MPI_Status*) calloc(mesh->nprocs, sizeof(MPI_Status));
    MPI_Waitall(Nmess, mpi_out_requests, outstatus);
    free(outstatus);

    free(mpi_out_requests);
    free(mpi_in_requests);
}