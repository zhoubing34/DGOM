#include <MultiRegions/MultiRegions.h>
#include <MultiRegions/MultiRegBC/MultiRegBC2d.h>
#include "VarMapPair_test.h"

void FetchParmapNode2d(MultiRegBC2d *mesh,
                       double *f_Q, double *f_inQ, double *f_outQ,
                       MPI_Request *mpi_send_requests,
                       MPI_Request *mpi_recv_requests, int *Nmess);

int MultiTriRegions_VarMapPair_Test(MultiRegBC2d *surf, int verbose){
    // regions
    int fail = 0, i;
//    MultiReg2d* mesh = setTriTestMesh();
    MultiReg2d *mesh = surf->mesh;
    StdRegions2d *shape = mesh->stdcell;

    // local variables
    int Nfp = mesh->K*shape->Nfp*shape->Nfaces;
    double *xM = Vector_create(Nfp);
    double *xP = Vector_create(Nfp);
    double *yM = Vector_create(Nfp);
    double *yP = Vector_create(Nfp);
    int parNodeTotalOut = surf->parNodeTotalOut;
    double *x_in  = Vector_create(parNodeTotalOut);
    double *x_out = Vector_create(parNodeTotalOut);
    double *y_in  = Vector_create(parNodeTotalOut);
    double *y_out = Vector_create(parNodeTotalOut);

    // local varmap
    int Nmess = 0;
    for(i=0; i<Nfp; i++) {
        int idM = surf->vmapM[i];

        xM[i] = mesh->x[0][idM];
        yM[i] = mesh->y[0][idM];
    }

    MPI_Request *mpi_out_requests = (MPI_Request*) calloc(mesh->nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_in_requests  = (MPI_Request*) calloc(mesh->nprocs, sizeof(MPI_Request));
    FetchParmapNode2d(surf, mesh->x[0], x_in, x_out, mpi_out_requests, mpi_in_requests, &Nmess);
    FetchParmapNode2d(surf, mesh->y[0], y_in, y_out, mpi_out_requests, mpi_in_requests, &Nmess);

    MPI_Status *instatus  = (MPI_Status*) calloc(mesh->nprocs, sizeof(MPI_Status));
    MPI_Waitall(Nmess, mpi_in_requests, instatus);
    free(instatus);

    for(i=0; i<Nfp; i++) {
        int idP = surf->vmapP[i];
        if (idP<0){
            int t = -idP-1;
            xP[i] = x_in[t];
            yP[i] = y_in[t];
        }else{
            xP[i] = mesh->x[0][idP];
            yP[i] = mesh->y[0][idP];
        }
    }

    if(!mesh->procid){
        fail = Vector_test("VarQuadMapPair_x", xP, xM, Nfp, 0);
        fail = Vector_test("VarQuadMapPair_y", yP, yM, Nfp, 0);
    }

    free(mpi_in_requests);
    free(mpi_out_requests);

    if(verbose){
        char casename[32] = "MultiTriRegions_VarMapPair_Test";
        FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);
        PrintVector2File(fp, "xM", xM, Nfp);
        PrintVector2File(fp, "xP", xP, Nfp);
        PrintVector2File(fp, "yM", yM, Nfp);
        PrintVector2File(fp, "yP", yP, Nfp);
        /* write vmapM */
        PrintIntVector2File(fp, "vmapM", surf->vmapM, Nfp);
        /* write vmapP */
        PrintIntVector2File(fp, "vmapP", surf->vmapP, Nfp);
        PrintIntVector2File(fp, "Npar", mesh->Npar, mesh->nprocs);
        fprintf(fp, "parNtotalout = %d\n", surf->parNodeTotalOut);
        PrintIntVector2File(fp, "parmapOut", surf->nodeIndexOut, surf->parNodeTotalOut);

        PrintVector2File(fp, "x_out", x_out, surf->parNodeTotalOut);
        PrintVector2File(fp, "y_out", y_out, surf->parNodeTotalOut);
        PrintVector2File(fp, "x_in", x_in, surf->parNodeTotalOut);
        PrintVector2File(fp, "y_in", y_in, surf->parNodeTotalOut);

        fclose(fp);
    }

    Vector_free(xM);
    Vector_free(xP);
    Vector_free(yM);
    Vector_free(yP);
    Vector_free(x_in);
    Vector_free(y_in);
    Vector_free(x_out);
    Vector_free(y_out);

    return fail;
}

int MultiQuadRegions_VarMapPair_Test(MultiRegBC2d *surf, int verbose){
    // local regions
    int i, fail = 0;
//    MultiReg2d* mesh = setQuadTestMesh();
    MultiReg2d *mesh = surf->mesh;
    StdRegions2d *shape = mesh->stdcell;

    // local variables
    int Nfp = mesh->K*shape->Nfp*shape->Nfaces;
    double *xM = Vector_create(Nfp);
    double *xP = Vector_create(Nfp);
    double *yM = Vector_create(Nfp);
    double *yP = Vector_create(Nfp);

    int parNodeTotalOut = surf->parNodeTotalOut;
    double *x_in  = Vector_create(parNodeTotalOut);
    double *x_out = Vector_create(parNodeTotalOut);
    double *y_in  = Vector_create(parNodeTotalOut);
    double *y_out = Vector_create(parNodeTotalOut);

    // local varmap
    int Nmess = 0;
    for(i=0; i<Nfp; i++) {
        int idM = surf->vmapM[i];

        xM[i] = mesh->x[0][idM];
        yM[i] = mesh->y[0][idM];
    }

    MPI_Request *mpi_out_requests = (MPI_Request*) calloc(mesh->nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_in_requests  = (MPI_Request*) calloc(mesh->nprocs, sizeof(MPI_Request));
    FetchParmapNode2d(surf, mesh->x[0], x_in, x_out, mpi_out_requests, mpi_in_requests, &Nmess);
    FetchParmapNode2d(surf, mesh->y[0], y_in, y_out, mpi_out_requests, mpi_in_requests, &Nmess);

    MPI_Status *instatus  = (MPI_Status*) calloc(mesh->nprocs, sizeof(MPI_Status));
    MPI_Waitall(Nmess, mpi_in_requests, instatus);
    free(instatus);

    for(i=0; i<Nfp; i++) {
        int idP = surf->vmapP[i];
        if (idP<0){
            int t = -idP-1;
            xP[i] = x_in[t];
            yP[i] = y_in[t];
        }else{
            xP[i] = mesh->x[0][idP];
            yP[i] = mesh->y[0][idP];
        }

    }

    if(!mesh->procid){
        fail = Vector_test("VarQuadMapPair_x", xP, xM, Nfp, 0);
        fail = Vector_test("VarQuadMapPair_y", yP, yM, Nfp, 0);
    }

    free(mpi_in_requests);
    free(mpi_out_requests);

    if(verbose){
        //check
        char casename[33] = "MultiQuadRegions_VarMapPair_Test";
        FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);
        PrintVector2File(fp, "xM", xM, Nfp);
        PrintVector2File(fp, "xP", xP, Nfp);
        PrintVector2File(fp, "yM", yM, Nfp);
        PrintVector2File(fp, "yP", yP, Nfp);
        /* write vmapM */
        PrintIntVector2File(fp, "vmapM", surf->vmapM, Nfp);
        /* write vmapP */
        PrintIntVector2File(fp, "vmapP", surf->vmapP, Nfp);
        PrintIntVector2File(fp, "Npar", mesh->Npar, mesh->nprocs);
        fprintf(fp, "parNtotalout = %d\n", surf->parNodeTotalOut);
        PrintIntVector2File(fp, "parmapOut", surf->nodeIndexOut,  surf->parNodeTotalOut);

        PrintVector2File(fp, "x_out", x_out, surf->parNodeTotalOut);
        PrintVector2File(fp, "y_out", y_out, surf->parNodeTotalOut);
        PrintVector2File(fp, "x_in", x_in, surf->parNodeTotalOut);
        PrintVector2File(fp, "y_in", y_in, surf->parNodeTotalOut);
        fclose(fp);
    }

    Vector_free(xM);
    Vector_free(xP);
    Vector_free(yM);
    Vector_free(yP);
    Vector_free(x_in);
    Vector_free(y_in);
    Vector_free(x_out);
    Vector_free(y_out);

    return fail;
}


void FetchParmapNode2d(MultiRegBC2d *surf,
                       double *f_Q, double *f_inQ, double *f_outQ,
                       MPI_Request *mpi_send_requests,
                       MPI_Request *mpi_recv_requests,
                       int *Nmessage) {

    MultiReg2d *mesh = surf->mesh;
    int t;
    /* buffer outgoing node data */
    for(t=0;t<surf->parNodeTotalOut;++t)
        f_outQ[t] = f_Q[surf->nodeIndexOut[t]];

    StdRegions2d *shape = mesh->stdcell;

    /* do sends */
    int sk = 0, Nmess = 0;
    int p, Nout;
    for(p=0;p<mesh->nprocs;++p){
        if(p!=mesh->procid){
            Nout = mesh->Npar[p]*shape->Nfp; // # of variables send to process p
            if(Nout){
                /* symmetric communications (different ordering) */
                MPI_Isend(f_outQ+sk, Nout, MPI_DOUBLE, p, 6666+p,
                          MPI_COMM_WORLD, mpi_send_requests +Nmess);
                MPI_Irecv(f_inQ+sk,  Nout, MPI_DOUBLE, p, 6666+mesh->procid,
                          MPI_COMM_WORLD,  mpi_recv_requests +Nmess);
                sk+=Nout;
                ++Nmess;
            }
        }
    }
    *Nmessage = Nmess; /* number of messages */
}