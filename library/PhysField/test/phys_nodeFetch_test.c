//
// Created by li12242 on 16/12/16.
//

#include <MultiRegions/mr_mesh.h>
#include "phys_nodeFetch_test.h"
#include "PhysField/phys_fetchBuffer.h"
#include "phys_test.h"


int phys_nodeFetch_test(physField *phys, int verbose, char *message, char *filename){
    // local variable
    int fail = 0;

    parallMesh *mesh = phys->mesh;
    multiReg *region = phys->region;
    stdCell *shape = phys->cell;

    const int K = phys->grid->K;
    const int Np = shape->Np;
    const int Nfp = shape->Nfp;
    const int Nfaces = shape->Nfaces;
    const int nprocs = mesh->nprocs;

    int Nnode = Nfaces*Nfp*K;
    double xM[Nnode];
    double xP[Nnode];
    double yM[Nnode];
    double yP[Nnode];

    int k,i;
    // assignment
    int sk = 0;
    for(k=0;k<K;k++){
        for(i=0;i<Np;i++){
            phys->f_Q[sk++] = region->x[k][i];
            phys->f_Q[sk++] = region->y[k][i];
        }
    }

    sk = 0;
    for(k=0;k<K;k++){
        for(i=0;i<Np*phys->Nfield;i++){
            phys->f_ext[sk] = phys->f_Q[sk];
            sk++;
        }
    }

    // MPI send & recv operations
    double clockT1, clockT2;
    int Nmess;
    MPI_Request mpi_send_requests[nprocs];
    MPI_Request mpi_recv_requests[nprocs];

//    printf("procid=%d, fetch node buffer\n", mesh->procid);
    clockT1 = MPI_Wtime();
    fetchNodeBuffer2d(phys, mpi_send_requests, mpi_recv_requests, &Nmess);
    clockT2 = MPI_Wtime();

    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_send_requests, instatus);
    MPI_Waitall(Nmess, mpi_recv_requests, instatus);

//    printf("procid=%d, finish MPI_Waitall\n", mesh->procid);

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

#if 0
            if(mesh->procid)
                printf("p=%d, k=%d, i=%d, bcType=%d, idM=%d, idP=%d\n",mesh->procid, k,i,bsType,idM,idP);
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

//    printf("procid=%d, finish xP,yP assignment\n", mesh->procid);

    if(!mesh->procid) {
        fail = Vector_test(message, xM, xP, Nnode, (clockT2 - clockT1));
        fail = Vector_test(message, yM, yP, Nnode, (clockT2 - clockT1));
    }

    if(verbose){
        FILE *fp = CreateLog(filename, mesh->procid, mesh->nprocs);
        fprintf(fp, "Nfield = %d\n", phys->Nfield);
        fprintf(fp, "dim = %d\n", phys->dim);
        fprintf(fp, "Nsurfinfo = %d\n", phys->Nsurfinfo);
        PrintVector2File(fp, "surinfo", phys->surfinfo, K*phys->Nsurfinfo*Nfp*Nfaces);
        fprintf(fp, "Nvgeo = %d\n", phys->Nvgeo);
        PrintVector2File(fp, "vgeo", phys->vgeo, K*phys->Nvgeo);
        PrintVector2File(fp, "f_Q", phys->f_Q, K*Np*phys->Nfield);
        PrintVector2File(fp, "xM", xM, Nnode);
        PrintVector2File(fp, "xP", xP, Nnode);
        PrintVector2File(fp, "yM", yM, Nnode);
        PrintVector2File(fp, "yP", yP, Nnode);
        fclose(fp);
    }

//    printf("step out phys_nodeFetch_test\n");

    return fail;
}