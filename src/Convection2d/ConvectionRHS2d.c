#include "ConvectionDriver2d.h"


void ConvectionRHS2d(PhysDomain2d *phys, PhysDomain2d *flowRate,
                     float frka, float frkb, float fdt){
    /* registers and temporary */
    register unsigned int k, t, geoid=0;

    MultiReg2d   *mesh  = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    /* mesh parameters */
    const int K = mesh->K;

    real *vgeo     = phys->vgeo;
    real *surfinfo = phys->surfinfo;
    real *f_Dr    = shape->f_Dr;
    real *f_Ds    = shape->f_Ds;
    real *f_LIFT  = shape->f_LIFT;

    real *f_Q     = phys->f_Q;
    real *f_rhsQ  = phys->f_rhsQ;
    real *f_resQ  = phys->f_resQ;
    real *f_s     = flowRate->f_Q; /* flow rate */

    real *f_inQ   = phys->f_inQ;

    /* mpi request buffer */
    MPI_Request *mpi_out_requests = (MPI_Request*) calloc(mesh->nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_in_requests  = (MPI_Request*) calloc(mesh->nprocs, sizeof(MPI_Request));

    int Nmess;

    FetchParmapNode2d(phys, mpi_out_requests, mpi_in_requests, &Nmess);

    /* NOTE: should be local memory */
    float  Qk[shape->Np*phys->Nfields];     /* scalar field */
    float  Uf[shape->Np*flowRate->Nfields]; /* flow rate field */

    /* volume integral */
    for(k=0;k<K;++k){

        /* NOTE: once k is known, all other indexing variables etc are derived */
        register unsigned int n, m;

        /* NOTE: buffer element k into local storage */
        float *qpt = f_Q + phys->Nfields*shape->Np*k;
        float *upt = f_s + flowRate->Nfields*shape->Np*k;

        int uk = 0;
        for(m=0;m<flowRate->Nfields*shape->Np;++m){
            Uf[uk++] = upt[uk];
        }

        for(m=0;m<phys->Nfields*shape->Np;++m){
            Qk[m] = qpt[m];
        }

        for(n=0;n<shape->Np;++n){

            const float *ptDr = f_Dr+n*shape->Np;
            const float *ptDs = f_Ds+n*shape->Np;

            const float drdx = vgeo[geoid++], drdy = vgeo[geoid++];
            const float dsdx = vgeo[geoid++], dsdy = vgeo[geoid++];

            float rhs = 0;

            int sk = 0; uk = 0;
            for(m=0;m<shape->Np;++m){
                const float dr = ptDr[m];
                const float ds = ptDs[m];
                const float dx = drdx*dr+dsdx*ds;
                const float dy = drdy*dr+dsdy*ds;
                const float u = Uf[uk++];
                const float v = Uf[uk++];

                const float C = Qk[sk++];

                rhs -= dx*u*C;
                rhs -= dy*v*C;
            }

            int id = phys->Nfields*(k*shape->Np + n);
            f_rhsQ[id] = rhs;
        }
    }

    /* DO RECV */
    MPI_Status *instatus  = (MPI_Status*) calloc(mesh->nprocs, sizeof(MPI_Status));
    MPI_Waitall(Nmess, mpi_in_requests, instatus);
    free(instatus);

    /* surface integral */
    for(k=0;k<K;++k){

        /* NOTE: once k is known, all other indexing variables etc are derived */
        register unsigned int n, m;

        /* NOTE: should be local memory */
        float fluxQ[shape->Nfaces*shape->Nfp*phys->Nfields];

        /* NOTE: index into geometric factors */
        int surfid=k*6*shape->Nfp*shape->Nfaces;

        /* Lax-Friedrichs flux */
        int sk = 0;
        for(m=0;m<shape->Nfp*shape->Nfaces;++m){

            int   idM       = (int)surfinfo[surfid++];
            int   idP       = (int)surfinfo[surfid++];
            const float FSc = surfinfo[surfid++];
            const float BSc = surfinfo[surfid++];
            const float NXf = surfinfo[surfid++];
            const float NYf = surfinfo[surfid++];

            float dC, dF, dG;
            if(idP<0){
                idP = phys->Nfields*(-1-idP);
                dC = (f_inQ[idP] - f_Q[idM]);
            }
            else if(idM == idP){
                dC =  0.0f - f_Q[idM];
            }else{
                dC = ( f_Q[idP] - f_Q[idM] );
            }

            const float uM = f_s[idM*flowRate->Nfields];
            const float vM = f_s[idM*flowRate->Nfields+1];

            dF = uM*dC;
            dG = vM*dC;

            const float un = fabsf(NXf*uM + NYf*vM);
            fluxQ[sk++] = FSc*( - NXf*dF - NYf*dG + un*dC)*(real)0.5;

        }

        /* Lift the flux term */
        for(n=0;n<shape->Np;++n){

            const float *ptLIFT = f_LIFT+n*shape->Nfp*shape->Nfaces;

            float rhs = 0;

            sk = 0;
            for(m=0;m<shape->Nfp*shape->Nfaces;++m){
                const float L = ptLIFT[m];
                rhs += L*fluxQ[sk++];
            }

            int id = phys->Nfields*(n+k*shape->Np);
            f_rhsQ[id] += rhs;
        }
    }


    for(t=0;t<K*shape->Np*phys->Nfields;++t){
        f_resQ[t] = frka*f_resQ[t]+fdt*f_rhsQ[t];   // calculate the resdiual of the equation
        f_Q[t]   += frkb*f_resQ[t];                 // evaluate scalar at next internal time step
    }

    /* make sure all messages went out */
    MPI_Status *outstatus  = (MPI_Status*) calloc(mesh->nprocs, sizeof(MPI_Status));
    MPI_Waitall(Nmess, mpi_out_requests, outstatus);
    free(outstatus);

    /* deallocate mem */
    free(mpi_out_requests);
    free(mpi_in_requests);
}