#include "SWEDriver2d.h"

void SWE_RHS2d(PhysDomain2d *phys, SWE_Solver2d *solver,
               const real frka, const real frkb, const real fdt){

    /* registers and temporary */
    register unsigned int k, geoid=0;

    MultiReg2d   *mesh  = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    /* mesh parameters */
    const int K       = mesh->K;
    const int Nfields = phys->Nfields;
    const int Np      = shape->Np;
    const int nprocs  = mesh->nprocs;

    real *vgeo     = phys->vgeo;
    real *surfinfo = phys->surfinfo;
    real *f_Dr    = shape->f_Dr;
    real *f_Ds    = shape->f_Ds;
    real *f_LIFT  = shape->f_LIFT;

    real *f_Q     = phys->f_Q;
    real *f_rhsQ  = phys->f_rhsQ;
    real *f_resQ  = phys->f_resQ;

    real *f_inQ   = phys->f_inQ;

    /* fetch boundary values */
    MPI_Request *mpi_out_requests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_in_requests  = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));

    int Nmess=0;
    FetchParmapNode2d(phys, mpi_out_requests, mpi_in_requests, &Nmess);

    /* volume integral */
    for(k=0;k<K;++k) {

        /* NOTE: once k is known, all other indexing variables etc are derived */
        register unsigned int n, m;

        /* NOTE: should be local memory */
        real Qk  [Np*Nfields];
        real Eflx[Np*Nfields];
        real Gflx[Np*Nfields];
        real Sour[Np*Nfields];

        /* NOTE: buffer element k into local storage */
        real *qpt = f_Q+Nfields*Np*k;
        for(m=0;m<Nfields*Np;++m){
            Qk[m] = qpt[m];
        }

        /* flux terms */
        SWE_ElementalFlux2d(phys, solver, Qk, Eflx, Gflx);
        SWE_ElementalSource2d(phys, solver, k, vgeo + k * Np * 4, Qk, Sour);

        for(n=0;n<Np;++n){

            const real *ptDr = f_Dr+n*Np;
            const real *ptDs = f_Ds+n*Np;

            const real drdx = vgeo[geoid++], drdy = vgeo[geoid++];
            const real dsdx = vgeo[geoid++], dsdy = vgeo[geoid++];

            real rhsH = 0, rhsQx = 0, rhsQy = 0;

            int sk = 0;
            for(m=0;m<Np;++m){
                const real dr = ptDr[m];
                const real ds = ptDs[m];
                const real dx = drdx*dr+dsdx*ds;
                const real dy = drdy*dr+dsdy*ds;

                rhsH  += -(dx*Eflx[sk]+dy*Gflx[sk]) + Sour[sk++];
                rhsQx += -(dx*Eflx[sk]+dy*Gflx[sk]) + Sour[sk++];
                rhsQy += -(dx*Eflx[sk]+dy*Gflx[sk]) + Sour[sk++];
            }

            int id = Nfields*(k*Np + n);
            f_rhsQ[id++] = rhsH;
            f_rhsQ[id++] = rhsQx;
            f_rhsQ[id]   = rhsQy;
        }
    }
    /* DO RECV */
    MPI_Status *instatus  = (MPI_Status*) calloc(nprocs, sizeof(MPI_Status));
    MPI_Waitall(Nmess, mpi_in_requests, instatus);
    free(instatus);


    const int Nfaces = shape->Nfaces;
    const int Nfp    = shape->Nfp;

    /* surface integral */
    for(k=0;k<K;++k){

        /* NOTE: once k is known, all other indexing variables etc are derived */
        register unsigned int n, m;
        /* NOTE: should be local memory */
        real fluxQ[Nfaces*Nfp*Nfields];

        /* NOTE: index into geometric factors */
        int surfid=k*6*Nfp*Nfaces;

        /* Lax-Friedrichs flux */
        int sk = 0;
        for(m=0;m<Nfp*Nfaces;++m){

            int   idM       = (int)surfinfo[surfid++];
            int   idP       = (int)surfinfo[surfid++];
            const real FSc = surfinfo[surfid++];
            const real BSc = surfinfo[surfid++];
            const real NXf = surfinfo[surfid++];
            const real NYf = surfinfo[surfid++];

            real hM,hP,qxM,qxP,qyM,qyP;
            real Fhs,Fqxs,Fqys;
            if(idP<0){
                idP = Nfields*(-1-idP);
                hM  = f_Q[idM++]; hP  = f_inQ[idP++];
                qxM = f_Q[idM++]; qxP = f_inQ[idP++];
                qyM = f_Q[idM  ]; qyP = f_inQ[idP  ];
            }else{
                hM  = f_Q[idM++]; hP  = f_Q[idP++];
                qxM = f_Q[idM++]; qxP = f_Q[idP++];
                qyM = f_Q[idM  ]; qyP = f_Q[idP  ];
            }
            real EhM,EqxM,EqyM;
            real GhM,GqxM,GqyM;

            SWE_NodalFlux2d(solver, hM, qxM, qyM, &EhM, &EqxM, &EqyM, &GhM, &GqxM, &GqyM);
            SWE_NodalNumFlux2d(solver, NXf, NYf, hM, hP, qxM, qxP, qyM, qyP, &Fhs, &Fqxs, &Fqys);

            Fhs  = EhM*NXf  + GhM*NYf  - Fhs;
            Fqxs = EqxM*NXf + GqxM*NYf - Fqxs;
            Fqys = EqyM*NXf + GqyM*NYf - Fqys;

            fluxQ[sk++] = FSc*Fhs;
            fluxQ[sk++] = FSc*Fqxs;
            fluxQ[sk++] = FSc*Fqys;

        }

        for(n=0;n<Np;++n){
            const real *ptLIFT = f_LIFT+n*Nfp*Nfaces;
            float rhsH = 0, rhsQx = 0, rhsQy = 0;
            sk = 0;
            for(m=0;m<Nfp*Nfaces;++m){
                const real L = ptLIFT[m];
                rhsH  += L*fluxQ[sk++];
                rhsQx += L*fluxQ[sk++];
                rhsQy += L*fluxQ[sk++];
            }

            int id = Nfields*(n+k*Np);
            f_rhsQ[id++] += rhsH;
            f_rhsQ[id++] += rhsQx;
            f_rhsQ[id]   += rhsQy;
        }
    }

    int n;
    for(n=0;n<K*Np*Nfields;++n){
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
