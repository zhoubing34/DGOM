//
// Created by li12242 on 12/11/16.
//

#include <MultiRegions/MultiRegions.h>
#include <MultiRegions/MultiRegBC/MultiRegBC2d.h>
#include "CellPair_test.h"

void FetchParmapEle2d(MultiRegBC2d *surf,
                      real *f_E, real *f_inE, real *f_outE,
                      MPI_Request *mpi_send_requests,
                      MPI_Request *mpi_recv_requests,
                      int *Nmessage);


int MultiTriRegions_CellPair_test(MultiRegBC2d *surf, int verbose){

    // local variable
    int i,k,fail = 0;

//    MultiReg2d* mesh = setTriTestMesh();
    MultiReg2d *mesh = surf->mesh;
    stdCell *shape = mesh->stdcell;

    double *xM = Vector_create(mesh->K);
    double *xP = Vector_create(mesh->K);
    double *yM = Vector_create(mesh->K);
    double *yP = Vector_create(mesh->K);
    int parCellTotalOut = surf->parCellTotalOut;
    double *x_in  = Vector_create(parCellTotalOut);
    double *x_out = Vector_create(parCellTotalOut);
    double *y_in  = Vector_create(parCellTotalOut);
    double *y_out = Vector_create(parCellTotalOut);

    // calculate cell centre
    for(k=0; k<mesh->K; k++){
        for(i=0; i<shape->Nfaces; i++){
            int idE = mesh->EToE[k][i];

            if(idE<0){
                int t1 = shape->Fmask[i][0] + k*shape->Np;
                int t2 = shape->Fmask[i][shape->Nfp-1] + k*shape->Np;

                xM[k] = (mesh->x[0][t1]+mesh->x[0][t2])/2; // face centre
                yM[k] = (mesh->y[0][t1]+mesh->y[0][t2])/2; // face centre
            }

        }
    }

    int Nmess = 0;
    MPI_Request *mpi_out_requests = (MPI_Request*) calloc(mesh->nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_in_requests  = (MPI_Request*) calloc(mesh->nprocs, sizeof(MPI_Request));
    FetchParmapEle2d(surf, xM, x_in, x_out, mpi_out_requests, mpi_in_requests, &Nmess);
    FetchParmapEle2d(surf, yM, y_in, y_out, mpi_out_requests, mpi_in_requests, &Nmess);

    MPI_Status *instatus  = (MPI_Status*) calloc(mesh->nprocs, sizeof(MPI_Status));
    MPI_Waitall(Nmess, mpi_in_requests, instatus);
    free(instatus);

    for(k=0; k<mesh->K; k++){
        for(i=0; i<shape->Nfaces; i++){
            int t = mesh->EToE[k][i];
            if(t<0){
                xP[k] = x_in[-t-1];
                yP[k] = y_in[-t-1];
            }
        }
    }

    if(!mesh->procid){
        fail = Vector_test("EleTriMapPair_x", xP, xM, parCellTotalOut, 0);
        fail = Vector_test("EleTriMapPair_y", yP, yM, parCellTotalOut, 0);
    }

    free(mpi_in_requests);
    free(mpi_out_requests);


    if(verbose){
        /* gen log filename */
        char casename[32] = "MultiTriRegions_CellPair_test";
        FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

        /* write EToV */
        PrintIntMatrix2File(fp, "EToV", mesh->EToV, mesh->K, shape->Nv);
        /* write EToE */
        PrintIntMatrix2File(fp, "EToE", mesh->EToE, mesh->K, shape->Nv);
        /* write EToF */
        PrintIntMatrix2File(fp, "EToF", mesh->EToF, mesh->K, shape->Nv);
        /* write EToP */
        PrintIntMatrix2File(fp, "EToP", mesh->EToP, mesh->K, shape->Nv);
        fprintf(fp, "parEtotalOut = %d\n", surf->parCellTotalOut);
        PrintIntVector2File(fp, "elemapOUT = ", surf->cellIndexOut, surf->parCellTotalOut);

        PrintVector2File(fp, "xM", xM, mesh->K);
        PrintVector2File(fp, "xP", xP, mesh->K);
        PrintVector2File(fp, "yM", yM, mesh->K);
        PrintVector2File(fp, "yP", yP, mesh->K);
        PrintVector2File(fp, "x_out", x_out, parCellTotalOut);
        PrintVector2File(fp, "y_out", y_out, parCellTotalOut);
        PrintVector2File(fp, "x_in", x_in, parCellTotalOut);
        PrintVector2File(fp, "y_in", y_in, parCellTotalOut);

        fclose(fp);
    }

    Vector_free(xM);
    Vector_free(yM);
    Vector_free(xP);
    Vector_free(yP);
    Vector_free(x_in);
    Vector_free(y_in);
    Vector_free(x_out);
    Vector_free(y_out);

    return fail;
}


int MultiQuadRegions_CellPair_test(MultiRegBC2d *surf, int verbose){
    // local variable
    int i,k,fail = 0;

//    MultiReg2d* mesh = setQuadTestMesh();
    MultiReg2d *mesh = surf->mesh;
    stdCell *shape = mesh->stdcell;

    double *xM = Vector_create(mesh->K);
    double *xP = Vector_create(mesh->K);
    double *yM = Vector_create(mesh->K);
    double *yP = Vector_create(mesh->K);
    int parCellTotalOut = surf->parCellTotalOut;
    double *x_in  = Vector_create(parCellTotalOut);
    double *x_out = Vector_create(parCellTotalOut);
    double *y_in  = Vector_create(parCellTotalOut);
    double *y_out = Vector_create(parCellTotalOut);

    // calculate cell centre
    for(k=0; k<mesh->K; k++){
        for(i=0; i<shape->Nfaces; i++){
            int idE = mesh->EToE[k][i];

            if(idE<0){
                int t1 = shape->Fmask[i][0] + k*shape->Np;
                int t2 = shape->Fmask[i][shape->Nfp-1] + k*shape->Np;

                xM[k] = (mesh->x[0][t1]+mesh->x[0][t2])/2; // face centre
                yM[k] = (mesh->y[0][t1]+mesh->y[0][t2])/2; // face centre
            }

        }
    }

    int Nmess = 0;
    MPI_Request *mpi_out_requests = (MPI_Request*) calloc(mesh->nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_in_requests  = (MPI_Request*) calloc(mesh->nprocs, sizeof(MPI_Request));
    FetchParmapEle2d(surf, xM, x_in, x_out, mpi_out_requests, mpi_in_requests, &Nmess);
    FetchParmapEle2d(surf, yM, y_in, y_out, mpi_out_requests, mpi_in_requests, &Nmess);

    MPI_Status *instatus  = (MPI_Status*) calloc(mesh->nprocs, sizeof(MPI_Status));
    MPI_Waitall(Nmess, mpi_in_requests, instatus);
    free(instatus);

    for(k=0; k<mesh->K; k++){
        for(i=0; i<shape->Nfaces; i++){
            int t = mesh->EToE[k][i];
            if(t<0){
                xP[k] = x_in[-t-1];
                yP[k] = y_in[-t-1];
            }
        }
    }

    if(!mesh->procid){
        fail = Vector_test("EleQuadMapPair_x", xP, xM, surf->parCellTotalOut, 0);
        fail = Vector_test("EleQuadMapPair_y", yP, yM, surf->parCellTotalOut, 0);
    }

    free(mpi_in_requests);
    free(mpi_out_requests);

    if(verbose){
        /* gen log filename */
        char casename[32] = "MultiQuadRegions_CellPair_test";
        FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);

        /* write EToV */
        PrintIntMatrix2File(fp, "EToV", mesh->EToV, mesh->K, shape->Nv);
        /* write EToE */
        PrintIntMatrix2File(fp, "EToE", mesh->EToE, mesh->K, shape->Nv);
        /* write EToF */
        PrintIntMatrix2File(fp, "EToF", mesh->EToF, mesh->K, shape->Nv);
        /* write EToP */
        PrintIntMatrix2File(fp, "EToP", mesh->EToP, mesh->K, shape->Nv);
        fprintf(fp, "parEtotalOut = %d\n", parCellTotalOut);
        PrintIntVector2File(fp, "elemapOUT = ", surf->cellIndexOut, parCellTotalOut);

        PrintVector2File(fp, "xM", xM, mesh->K);
        PrintVector2File(fp, "xP", xP, mesh->K);
        PrintVector2File(fp, "yM", yM, mesh->K);
        PrintVector2File(fp, "yP", yP, mesh->K);
        PrintVector2File(fp, "x_out", x_out, parCellTotalOut);
        PrintVector2File(fp, "y_out", y_out, parCellTotalOut);
        PrintVector2File(fp, "x_in", x_in, parCellTotalOut);
        PrintVector2File(fp, "y_in", y_in, parCellTotalOut);

        fclose(fp);
    }

    Vector_free(xM);
    Vector_free(yM);
    Vector_free(xP);
    Vector_free(yP);
    Vector_free(x_in);
    Vector_free(y_in);
    Vector_free(x_out);
    Vector_free(y_out);

    return fail;
}


void FetchParmapEle2d(MultiRegBC2d *surf,
                      real *f_E, real *f_inE, real *f_outE,
                      MPI_Request *mpi_send_requests,
                      MPI_Request *mpi_recv_requests,
                      int *Nmessage){

    /* buffer outgoing node data */
    MultiReg2d *mesh = surf->mesh;
    int n;
    for(n=0;n<surf->parCellTotalOut;++n)
        f_outE[n] = f_E[surf->cellIndexOut[n]];

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