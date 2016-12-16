//
// Created by li12242 on 16/12/16.
//

#include <MultiRegions/MultiRegions.h>
#include "NodeFetch_test.h"
#include "PhysDomain_test.h"

#define DEBUG 0

int nodeFetch_test(PhysDomain2d *phys, int verbose, char *message, char *filename){
    // local variable
    int fail = 0;

    MultiReg2d *mesh = phys->mesh;
    MultiRegBC2d *surf = phys->surf;
    stdCell *shape = mesh->stdcell;

    int K = mesh->K;
    int Np = shape->Np;
    int Nfp = shape->Nfp;
    int Nfaces = shape->Nfaces;
    int nprocs = mesh->nprocs;

    int Nnode = Nfaces*Nfp*K;
    double *xM = Vector_create(Nnode);
    double *xP = Vector_create(Nnode);
    double *yM = Vector_create(Nnode);
    double *yP = Vector_create(Nnode);

    int k,i;
    // assignment
    int sk = 0;
    for(k=0;k<K;k++){
        for(i=0;i<Np;i++){
            phys->f_Q[sk++] = mesh->x[k][i];
            phys->f_Q[sk++] = mesh->y[k][i];
        }
    }

    sk = 0;
    for(k=0;k<K;k++){
        for(i=0;i<Np*phys->Nfields;i++){
            phys->f_ext[sk] = phys->f_Q[sk];
            sk++;
        }
    }

    // MPI send & recv operations
    clock_t clockT1, clockT2;
    int Nmess = 0;
    MPI_Request *mpi_send_requests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_recv_requests  = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));

    clockT1 = clock();
    fetchNodeBuffer2d(phys, mpi_send_requests, mpi_recv_requests, &Nmess);
    clockT2 = clock();

    MPI_Status *instatus  = (MPI_Status*) calloc(nprocs, sizeof(MPI_Status));
    MPI_Waitall(Nmess, mpi_recv_requests, instatus);
    free(instatus);

    sk = 0;
    for(k=0;k<K;k++){
        int surfid = k*phys->Nsurfinfo*Nfp*Nfaces;
        for(i=0;i<Nfaces*Nfp;i++){
            int idM = (int)phys->surfinfo[surfid++];
            int idP = (int)phys->surfinfo[surfid++];
            surfid++;
            int bsType = (int)phys->surfinfo[surfid++];
            surfid++;
            surfid++;

#if DEBUG
            if(!mesh->procid)
                printf("k=%d, i=%d, bcType=%d, idM=%d, idP=%d\n",k,i,bsType,idM,idP);
#endif

            xM[sk] = phys->f_Q[idM++];
            yM[sk] = phys->f_Q[idM++];

            switch (bsType){
                case INNERLOC:
                    xP[sk] = phys->f_Q[idP++];
                    yP[sk] = phys->f_Q[idP++];
                    break;
                case INNERBS:
                    xP[sk] = phys->f_inQ[idP++];
                    yP[sk] = phys->f_inQ[idP++];
                    break;
                default: // open boundary
                    xP[sk] = phys->f_ext[idP++];
                    yP[sk] = phys->f_ext[idP++];
                    break;
            }
            sk++;
        }
    }

    if(!mesh->procid) {
        fail = Vector_test(message, xM, xP, Nnode, (clockT2 - clockT1) / CLOCKS_PER_SEC);
        fail = Vector_test(message, yM, yP, Nnode, (clockT2 - clockT1) / CLOCKS_PER_SEC);
    }

    if(verbose){
        FILE *fp = CreateLog(filename, mesh->procid, mesh->nprocs);
        PrintVector2File(fp, "surinfo", phys->surfinfo, K*phys->Nsurfinfo*Nfp*Nfaces);
        PrintVector2File(fp, "f_Q", phys->f_Q, K*Np*phys->Nfields);
        PrintVector2File(fp, "xM", xM, Nnode);
        PrintVector2File(fp, "xP", xP, Nnode);
        PrintVector2File(fp, "yM", yM, Nnode);
        PrintVector2File(fp, "yP", yP, Nnode);
        fclose(fp);
    }

    Vector_free(xM);
    Vector_free(xP);
    Vector_free(yM);
    Vector_free(yP);

    free(mpi_send_requests);
    free(mpi_recv_requests);

    return fail;
}