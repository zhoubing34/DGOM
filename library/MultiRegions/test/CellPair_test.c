//
// Created by li12242 on 12/11/16.
//

#include "CellPair_test.h"

void FetchParmapEle2d(MultiReg2d *mesh,
                      real *f_E, real *f_inE, real *f_outE,
                      MPI_Request *mpi_send_requests,
                      MPI_Request *mpi_recv_requests,
                      int *Nmessage);


int MultiTriRegions_CellPair_test(MultiReg2d *mesh, int verbose){

    // local variable
    int i,k,fail = 0;

//    MultiReg2d* mesh = SetTriParallelMultiRegions();
    StdRegions2d *shape = mesh->stdcell;

    double *xM = Vector_create(mesh->K);
    double *xP = Vector_create(mesh->K);
    double *yM = Vector_create(mesh->K);
    double *yP = Vector_create(mesh->K);
    double *x_in  = Vector_create(mesh->parEtotalout);
    double *x_out = Vector_create(mesh->parEtotalout);
    double *y_in  = Vector_create(mesh->parEtotalout);
    double *y_out = Vector_create(mesh->parEtotalout);

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
    FetchParmapEle2d(mesh, xM, x_in, x_out, mpi_out_requests, mpi_in_requests, &Nmess);
    FetchParmapEle2d(mesh, yM, y_in, y_out, mpi_out_requests, mpi_in_requests, &Nmess);

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
        fail = Vector_test("EleTriMapPair_x", xP, xM, mesh->parEtotalout, 0);
        fail = Vector_test("EleTriMapPair_y", yP, yM, mesh->parEtotalout, 0);
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
        fprintf(fp, "parEtotalOut = %d\n", mesh->parEtotalout);
        PrintIntVector2File(fp, "elemapOUT = ", mesh->elemapOut, mesh->parEtotalout);

        PrintVector2File(fp, "xM", xM, mesh->K);
        PrintVector2File(fp, "xP", xP, mesh->K);
        PrintVector2File(fp, "yM", yM, mesh->K);
        PrintVector2File(fp, "yP", yP, mesh->K);
        PrintVector2File(fp, "x_out", x_out, mesh->parEtotalout);
        PrintVector2File(fp, "y_out", y_out, mesh->parEtotalout);
        PrintVector2File(fp, "x_in", x_in, mesh->parEtotalout);
        PrintVector2File(fp, "y_in", y_in, mesh->parEtotalout);

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

    // free and finalize
//    StdRegions2d_free(shape);
//    MultiReg2d_free(mesh);
    return fail;
}


int MultiQuadRegions_CellPair_test(MultiReg2d *mesh, int verbose){
    // local variable
    int i,k,fail = 0;

//    MultiReg2d* mesh = SetQuadParallelMultiRegions();
    StdRegions2d *shape = mesh->stdcell;

    double *xM = Vector_create(mesh->K);
    double *xP = Vector_create(mesh->K);
    double *yM = Vector_create(mesh->K);
    double *yP = Vector_create(mesh->K);
    double *x_in  = Vector_create(mesh->parEtotalout);
    double *x_out = Vector_create(mesh->parEtotalout);
    double *y_in  = Vector_create(mesh->parEtotalout);
    double *y_out = Vector_create(mesh->parEtotalout);

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
    FetchParmapEle2d(mesh, xM, x_in, x_out, mpi_out_requests, mpi_in_requests, &Nmess);
    FetchParmapEle2d(mesh, yM, y_in, y_out, mpi_out_requests, mpi_in_requests, &Nmess);

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
        fail = Vector_test("EleQuadMapPair_x", xP, xM, mesh->parEtotalout, 0);
        fail = Vector_test("EleQuadMapPair_y", yP, yM, mesh->parEtotalout, 0);
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
        fprintf(fp, "parEtotalOut = %d\n", mesh->parEtotalout);
        PrintIntVector2File(fp, "elemapOUT = ", mesh->elemapOut, mesh->parEtotalout);

        PrintVector2File(fp, "xM", xM, mesh->K);
        PrintVector2File(fp, "xP", xP, mesh->K);
        PrintVector2File(fp, "yM", yM, mesh->K);
        PrintVector2File(fp, "yP", yP, mesh->K);
        PrintVector2File(fp, "x_out", x_out, mesh->parEtotalout);
        PrintVector2File(fp, "y_out", y_out, mesh->parEtotalout);
        PrintVector2File(fp, "x_in", x_in, mesh->parEtotalout);
        PrintVector2File(fp, "y_in", y_in, mesh->parEtotalout);

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

    // free and finalize
//    StdRegions2d_free(shape);
//    MultiReg2d_free(mesh);
    return fail;
}


void FetchParmapEle2d(MultiReg2d *mesh,
                      real *f_E, real *f_inE, real *f_outE,
                      MPI_Request *mpi_send_requests,
                      MPI_Request *mpi_recv_requests,
                      int *Nmessage){

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