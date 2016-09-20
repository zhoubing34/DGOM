//
// Created by li12242 on 16/8/2.
//

#include <PhysDomain/PhysDomain.h>
#include "UnitTest/MultiRegions/MultiRegionsTest.h"
#include "SWE2d/SWEDriver2d.h"

/* Private function */
void PrintVec(char *message, real *vec, int n);
void PrintVar(PhysDomain2d *phys);
void FluxTermTest(PhysDomain2d *phys, SWE_Solver2d *solver);
void NumFLuxTest(PhysDomain2d *phys, SWE_Solver2d *solver);
void VolIntegralTest(PhysDomain2d *phys, SWE_Solver2d *solver);
void SurfIntegralTest(PhysDomain2d *phys, SWE_Solver2d *solver);

int main(int argc, char **argv){

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    /* set mesh and physical domain */
    int N=1, Nfields=3;
    StdRegions2d *tri = GenStdTriEle(N);
    MultiReg2d   *triMesh;
    SetTestTriMesh(tri, triMesh);
    PhysDomain2d *phys = GenPhysDomain2d(triMesh, Nfields);
    int i,k;
    for(k=0;k<triMesh->K;k++){
        for(i=0;i<tri->Np;i++) {
            int ind = (k*tri->Np + i)*Nfields;
            /* initial value assignment */
            if(k<=1){
                phys->f_Q[ind++] = 10.0;
                phys->f_Q[ind++] = 2.0;
                phys->f_Q[ind++] = 0;
            }else{
                phys->f_Q[ind++] = 2.0;
                phys->f_Q[ind++] = 1.0;
                phys->f_Q[ind++] = 0;
            }
        }
    }

    printf("Initial conditions\n");
    PrintVar(phys);

    /* solver */
    SWE_Solver2d *solver = (SWE_Solver2d*) calloc(1, sizeof(SWE_Solver2d));
    solver->hcrit = 1.0e-4;
    solver->gra   = 9.81;
    solver->bot   = BuildMatrix(triMesh->K, tri->Np); /* set bottom topography */

    for(k=0;k<triMesh->K;k++){
        for(i=0;i<tri->Np;i++)
            solver->bot[k][i] = triMesh->x[k][i];
    }

    FluxTermTest(phys, solver);
    NumFLuxTest(phys, solver);
//    VolIntegralTest(phys, solver);
//    SurfIntegralTest(phys, solver);

    MPI_Finalize();

    return 0;
}

void SurfIntegralTest(PhysDomain2d *phys, SWE_Solver2d *solver){

    register unsigned int k;

    MultiReg2d   *mesh  = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    /* mesh parameters */
    const int K       = mesh->K;
    const int Nfields = phys->Nfields;
    const int Np      = shape->Np;
    const size_t nprocs  = (size_t)mesh->nprocs;

    real *vgeo     = phys->vgeo;
    real *surfinfo = phys->surfinfo;
    real *f_Dr    = shape->f_Dr;
    real *f_Ds    = shape->f_Ds;
    real *f_LIFT  = shape->f_LIFT;

    real *f_Q     = phys->f_Q;
    real *f_rhsQ  = phys->f_rhsQ;
    real *f_resQ  = phys->f_resQ;

    real *f_inQ   = phys->f_inQ;
    real *f_outQ  = phys->f_outQ;

    const int Nfaces = shape->Nfaces;
    const int Nfp    = shape->Nfp;

    /* fetch boundary values */
    MPI_Request *mpi_out_requests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_in_requests  = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));

    int Nmess;
    FetchParmapNode2d(phys, mpi_out_requests, mpi_in_requests, &Nmess);
    printf("procid=%d, finish fetch nodes values\n",mesh->procid);

    /* DO RECV */
    MPI_Status *instatus  = (MPI_Status*) calloc(nprocs, sizeof(MPI_Status));
    MPI_Waitall(Nmess, mpi_in_requests, instatus);
    free(instatus);

    /* surface integral */
    for(k=0;k<K;++k){

        /* NOTE: once k is known, all other indexing variables etc are derived */
        register unsigned int n, m;
        /* NOTE: should be local memory */
        float fluxQ[Nfaces*Nfp*Nfields];

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

            printf("dFh=%f, dFqx=%f, dFqy=%f\n", Fhs, Fqxs, Fqys);

            fluxQ[sk++] = FSc*Fhs;
            fluxQ[sk++] = FSc*Fqxs;
            fluxQ[sk++] = FSc*Fqys;
        }

        for(n=0;n<Np;++n){
            const float *ptLIFT = f_LIFT+n*Nfp*Nfaces;
            float rhsH = 0, rhsQx = 0, rhsQy = 0;
            sk = 0;
            for(m=0;m<Nfp*Nfaces;++m){
                const float L = ptLIFT[m];
                rhsH +=  L*fluxQ[sk++];
                rhsQx += L*fluxQ[sk++];
                rhsQy += L*fluxQ[sk++];
            }

            int id = Nfields*(n+k*Np);
            f_rhsQ[id++] += rhsH;
            f_rhsQ[id++] += rhsQx;
            f_rhsQ[id]   += rhsQy;
        }
    }

    /* make sure all messages went out */
    MPI_Status *outstatus  = (MPI_Status*) calloc(mesh->nprocs, sizeof(MPI_Status));
    MPI_Waitall(Nmess, mpi_out_requests, outstatus);
    free(outstatus);

    free(mpi_out_requests);
    free(mpi_in_requests);
}

void VolIntegralTest(PhysDomain2d *phys, SWE_Solver2d *solver){

    /* registers and temporary */
    register unsigned int k, geoid=0;

    MultiReg2d   *mesh  = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    /* mesh parameters */
    const int K       = mesh->K;
    const int Nfields = phys->Nfields;
    const int Np      = shape->Np;

    real *vgeo    = phys->vgeo;
    real *f_Dr    = shape->f_Dr;
    real *f_Ds    = shape->f_Ds;
    real *f_Q     = phys->f_Q;
    real *f_rhsQ  = phys->f_rhsQ;

    printf("Volume integral test\n");

    for(k=0;k<K;k++){
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

            printf("drdx=%f,drdy=%f,dsdx=%f,dsdy=%f\n",drdx,drdy,dsdx,dsdy);

            real rhsH = 0, rhsQx = 0, rhsQy = 0;

            int sk = 0;
            for(m=0;m<Np;++m){
                const real dr = ptDr[m];
                const real ds = ptDs[m];
                const real dx = drdx*dr+dsdx*ds;
                const real dy = drdy*dr+dsdy*ds;

                rhsH  += +dx*Eflx[sk]+dy*Gflx[sk] + Sour[sk++];
                rhsQx += +dx*Eflx[sk]+dy*Gflx[sk] + Sour[sk++];
                rhsQy += +dx*Eflx[sk]+dy*Gflx[sk] + Sour[sk++];
            }

            int id = Nfields*(k*Np + n);
            f_rhsQ[id++] = rhsH;
            f_rhsQ[id++] = rhsQx;
            f_rhsQ[id]   = rhsQy;

            printf("k=%d, n=%d, rhsH=%f, rhsQx=%f, rhsQy=%f\n",
                   k,n,rhsH,rhsQx,rhsQy);
        }

        int sk = 0;
        for(n=0;n<Np;n++){
            printf("n =%2d, Eh =%f, Eqx =%f, Eqy =%f, Gh =%f, Gqx =%f, Gqy =%f, Sh=%f, Sqx=%f, Sqy=%f\n",
                   n, Eflx[sk], Eflx[sk+1], Eflx[sk+2], Gflx[sk], Gflx[sk+1], Gflx[sk+2], Sour[sk], Sour[sk+1], Sour[sk+2]); sk +=3;
        }
    }
}

void PrintVar(PhysDomain2d *phys){
    MultiReg2d   *mesh  = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    real *var = (real *) malloc(mesh->K*shape->Np*sizeof(real));
    int k,i,ind,sk=0;
    int Np = shape->Np, Nfields=phys->Nfields;
    for (k=0;k<mesh->K;k++){
        for (i=0;i<Np;i++){
            ind = (k*Np + i)*Nfields;
            // printf("k = %d, i = %d, ind = %d\n", k, i, ind);
            var[sk++] = phys->f_Q[ind];
        }
        PrintVec("h ",var+k*Np, Np);
    }

    sk = 0;
    for (k=0;k<mesh->K;k++){
        for (i=0;i<Np;i++){
            ind = (k*Np + i)*Nfields + 1;
            // printf("k = %d, i = %d, ind = %d\n", k, i, ind);
            var[sk++] = phys->f_Q[ind];
        }
        PrintVec("qx",var+k*Np, Np);
    }

    sk = 0;
    for (k=0;k<mesh->K;k++){
        for (i=0;i<Np;i++){
            ind = (k*Np + i)*Nfields + 2;
            // printf("k = %d, i = %d, ind = %d\n", k, i, ind);
            var[sk++] = phys->f_Q[ind];
        }
        PrintVec("qy",var+k*Np, Np);
    }

    free(var);
}

void PrintVec(char *message, real *vec, int n){
    int i;
    printf("%s: ", message);
    for (i=0;i<n;i++){
        printf("%f, ", vec[i]);
    }
    printf("\n");
}

void FluxTermTest(PhysDomain2d *phys, SWE_Solver2d *solver){

    MultiReg2d   *mesh  = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    int Nfields = phys->Nfields;
    int Np      = shape->Np;

    real Qk  [Np*Nfields];
    real Eflx[Np*Nfields];
    real Gflx[Np*Nfields];

    int k, m;

    for (k=0;k<mesh->K;k++) {
        real *qpt = phys->f_Q + Nfields * Np * k;
        for (m = 0; m < Nfields * Np; ++m) {
            Qk[m] = qpt[m];
        }

        SWE_ElementalFlux2d(phys, solver, Qk, Eflx, Gflx);
        int sk=0;
        for(m=0;m<Np;m++){
            printf("Eh=%f, Eqx=%f, Eqy=%f, Gh=%f, Gqx=%f, Gqy=%f\n",
                   Eflx[sk],Eflx[sk+1],Eflx[sk+2],Gflx[sk],Gflx[sk+1],Gflx[sk+2]);
            sk += 3;
        }
    }
}

void NumFLuxTest(PhysDomain2d *phys, SWE_Solver2d *solver){

    MultiReg2d   *mesh  = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    const int Nfaces = shape->Nfaces;
    const int Nfp    = shape->Nfp;
    const int K      = mesh->K;
    int   Nfields  = phys->Nfields;
    real *surfinfo = phys->surfinfo;
    real *f_Q      = phys->f_Q;
    real *f_inQ    = phys->f_inQ;

    int k;
    /* surface integral */
    for(k=0;k<K;++k){

        /* NOTE: once k is known, all other indexing variables etc are derived */
        register unsigned int n, m;
        /* NOTE: should be local memory */
        float fluxQ[Nfaces*Nfp*Nfields];

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
                qyM = f_Q[idM++]; qyP = f_inQ[idP++];
            }else{
                hM  = f_Q[idM++]; hP  = f_Q[idP++];
                qxM = f_Q[idM++]; qxP = f_Q[idP++];
                qyM = f_Q[idM++]; qyP = f_Q[idP++];
            }

//            printf("ele:%d, Nfp=%d, ", k, m);
            SWE_NodalNumFlux2d(solver, NXf, NYf,
                               hM, hP, qxM, qxP, qyM, qyP,
                               &Fhs, &Fqxs, &Fqys);

            fluxQ[sk++] = FSc*Fhs;
            fluxQ[sk++] = FSc*Fqxs;
            fluxQ[sk++] = FSc*Fqys;

            printf("ele:%d, Nfp:%d, Fhs:%f, Fqxs:%f, Fqy:%f\n", k, m, Fhs, Fqxs, Fqys);
        }
    }
}