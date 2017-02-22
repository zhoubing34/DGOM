//
// Created by li12242 on 16/12/16.
//

#include <MultiRegions/mr_mesh.h>
#include "pf_nodeFetch_test.h"
#include "PhysField/pf_fetchBuffer.h"
#include "pf_test.h"

int phys_nodeFetch_test(physField *phys, int verbose){
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
            phys->f_Q[sk] = region->x[k][i];
            phys->f_ext[sk] = phys->f_Q[sk];
            sk++;

            phys->f_Q[sk] = region->y[k][i];
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
    pf_fetchNodeBuffer2d(phys, mpi_send_requests, mpi_recv_requests, &Nmess);
    clockT2 = MPI_Wtime();

    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_send_requests, instatus);
    MPI_Waitall(Nmess, mpi_recv_requests, instatus);

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
        fail = Vector_test(__FUNCTION__, xM, xP, Nnode, (clockT2 - clockT1));
        fail = Vector_test(__FUNCTION__, yM, yP, Nnode, (clockT2 - clockT1));
    }

    if(verbose){
        FILE *fp = CreateLog(__FUNCTION__, mesh->procid, mesh->nprocs);
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

    return fail;
}