//
// Created by li12242 on 16/8/3.
//

#include <MultiRegions/mr_mesh.h>
#include "phys_limiter.h"
#include "phys_cellMean.h"
#include "phys_fetchBuffer.h"

#define DEBUG 0
#if DEBUG
#include "LibUtilities/UTest.h"
#endif

#define TOL 1e-8
#define EXCEPTION(u, u1, u2) ( ((u-u1)>TOL)& ((u-u2)>TOL) )||( ((u1-u)>TOL)&((u2-u)>TOL))

/**
 * @brief calculate the gradient by Green method
 *
 * @param[in] phys
 * @param[in] uf_mean
 * @param[out] px
 * @param[out] py
 */
static void phys_gradient_Green(physField *phys, real *uf_mean, real *px, real *py){
    const int K = phys->grid->K;
    const int Nfield = phys->Nfield;
    const int Nfaces = phys->cell->Nfaces;
    const int Nfp = phys->cell->Nfp;

    double **x = phys->region->x;
    double **y = phys->region->y;

    register int k,f,fld,sk=0;
    for(k=0;k>K;k++){
        double A = 1.0/phys->region->size[k];
        for(fld=0;fld<Nfield;fld++){
            px[sk] = 0.;
            py[sk] = 0.;

            for(f=0;f<Nfaces;f++){ // loop over each faces
                int v1 = phys->cell->Fmask[f][0];
                int v2 = phys->cell->Fmask[f][Nfp-1];
                double x1 = x[k][v1];
                double x2 = x[k][v2];
                double y1 = y[k][v1];
                double y2 = y[k][v2];

                double dx =  (y2-y1);
                double dy = -(x2-x1);

                real uf = uf_mean[(k*Nfaces+f)*Nfield+fld];
                px[sk] += uf*dx*A;
                py[sk] += uf*dy*A;
            }
        }

    }
    return;
}

/**
 * @brief calculate the averaged mean value of each face
 * @param[in] phys
 * @param[out] uf_mean
 */
static void phys_faceMean(physField *phys, real *uf_mean){
    const int K = phys->grid->K;
    const int Nfield = phys->Nfield;
    const int Nfaces = phys->cell->Nfaces;
    const int Nfp = phys->cell->Nfp;
    const int Np = phys->cell->Np;

    double *ws = phys->cell->ws;

    register int k,f,fld,n;
    int sk=0;
    for(k=0;k<K;k++){
        real *f_Q = phys->f_Q + k*Np*Nfield; // values in the k-th cell
        for(f=0;f<Nfaces;f++){
            int *fmask = phys->cell->Fmask[f]; // face node index

            for(fld=0;fld<Nfield;fld++){
                uf_mean[sk] = 0.;
                for(n=0;n<Nfp;n++){
                    uf_mean[sk] += f_Q[ fmask[n]*Nfield+fld ]*ws[n]*0.5;
                }
                sk++;
            }
        }
    }
    return;
}

/**
 * @brief
 * Slope limiter from Anastasiou and Chan (1997) for two dimensional problems.
 *
 * @details
 * The slope limiter will act on each physical variables and reconstruct the
 * scalar distribution.
 *
 * @param[in] phys
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * phys   |   PhysDomain2d*  | contains reconstruct values
 *
 */
void phys_slloc2d(physField *phys, double beta){

    const int K = phys->grid->K;
    const int Nfield = phys->Nfield;
    const int procid = phys->mesh->procid;
    const int nprocs = phys->mesh->nprocs;
    const int Nfaces = phys->cell->Nfaces;
    const int Np = phys->cell->Np;

    register int k,n,f,fld;

#if DEBUG
    FILE *fp = CreateLog("phys_slloc2d", phys->mesh->procid, phys->mesh->nprocs);
#endif

    /* 1. calculate the cell average value */
    phys_cellMean(phys);

    /* 2. fetch cell info with other processes */
    MPI_Request mpi_out_requests[nprocs];
    MPI_Request mpi_in_requests[nprocs];
    int Nmess;

    /* do sends and recv */
    phys_fetchCellBuffer(phys, mpi_out_requests, mpi_in_requests, &Nmess);

    /* 3. face averaged values */
    real uf_mean[K*Nfaces*Nfield];
    phys_faceMean(phys, uf_mean);

#if DEBUG
    PrintVector2File(fp, "uf_mean", uf_mean, K*Nfaces*Nfield);
#endif

    /* 4.1. judge the trouble cell, get c_max and c_min  */
    int **EToE = phys->mesh->EToE;
    int **EToP = phys->mesh->EToP;
    real cell_max[K*Nfield], cell_min[K*Nfield];
    int trouble_cell[K],sk=0;

    for(k=0;k<K;k++){
        trouble_cell[k] = 0;
        real *c_mean = phys->c_Q + k*Nfield; // mean value of k-th cell
        for(fld=0;fld<Nfield;fld++){
            real temp = c_mean[fld]; // mean value of fld physical field
            /* initialize with k-th mean value */
            real c_max = temp;
            real c_min = temp;

            /* loop over all faces */
            for(f=0;f<Nfaces;f++){
                int e = EToE[k][f];
                int p = EToP[k][f];

                // if adjacent cell is on other process, jump this cycle
                if(p!=procid) continue;
                real c_next = phys->c_Q[e*Nfield+fld];

                c_max = max(c_max, c_next);
                c_min = min(c_min, c_next);

                if( EXCEPTION(uf_mean[(k*Nfaces+f)*Nfield+fld], temp, c_next) ){
#if DEBUG
                    fprintf(fp, "k=%d, f=%d, fld=%d, c_mean=%f, c_next=%f, uf_mean=%f\n",
                            k, f, fld, temp, c_next, uf_mean[(k*Nfaces+f)*Nfield+fld]);
#endif
                    trouble_cell[k] = 1;
                }
            }
            cell_max[sk  ] = c_max;
            cell_min[sk++] = c_min;
        }
    }

    /* wait for 2. sends and recv buffers */
    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_in_requests, instatus);
    MPI_Waitall(Nmess, mpi_out_requests, instatus);

    /* 4.2. judge the trouble cell, get c_max and c_min from parallel process  */
    parallMesh *mesh = phys->mesh;
    for(n=0;n<mesh->parallCellNum;n++){
        k = mesh->cellIndexIn[n];
        f = mesh->faceIndexIn[n];
        real *c_mean = phys->c_Q+k*Nfield; // mean value of k-th cell
        for(fld=0;fld<Nfield;fld++){
            real temp = c_mean[fld];
            real c_next = phys->c_inQ[n*Nfield+fld];

            cell_max[k*Nfield+fld] = max(cell_max[k*Nfield+fld], c_next);
            cell_min[k*Nfield+fld] = min(cell_min[k*Nfield+fld], c_next);
            if( EXCEPTION(uf_mean[(k*Nfaces+f)*Nfield+fld], temp, c_next) ){
#if DEBUG
                fprintf(fp, "k=%d, f=%d, fld=%d, c_mean=%f, c_next=%f, uf_mean=%f\n",
                        k, f, fld, temp, c_next, uf_mean[(k*Nfaces+f)*Nfield+fld]);
#endif
                trouble_cell[k] = 1;
            }
        }
    }

#if DEBUG
    PrintIntVector2File(fp, "trouble_cell", trouble_cell, K);
#endif

    /* 5. calculate the unlimited gradient */
    real px[K*Nfield], py[K*Nfield];
    phys_gradient_Green(phys, uf_mean, px, py);

    /* 6. calculate the limited results */
    multiReg *region = phys->region;
    for(k=0;k<K;k++){
        /* check if the cell need limited */
        if(trouble_cell[k]==0 ) continue;

        double xc = mr_reg_integral(region, k, region->x[k]);
        double yc = mr_reg_integral(region, k, region->y[k]);

        real *f_Q = phys->f_Q + k*Np*Nfield; // variable of k-th cell
        for(fld=0;fld<Nfield;fld++){
            real qmax = cell_max[k*Nfield+fld];
            real qmin = cell_min[k*Nfield+fld];
            real qmean = phys->c_Q[k*Nfield+fld];
            real t, psi = 1.0;
            for(n=0;n<Np;n++){
                real qval = f_Q[n*Nfield+fld];
                double rk = 1.0;
                /* compute the limiter `psi` */
                if(qval > qmax){
                    rk = (qmax - qmean)/(qval - qmean);
                }else if(qval < qmin){
                    rk = (qmin - qmean)/(qval - qmean);
                }
                t   = max(min(beta*rk, 1.0), min(rk, beta));
                psi = min(psi, t);
            }

            /* compute the limited gradient */
            real qx = px[k*Nfield+fld] * psi;
            real qy = py[k*Nfield+fld] * psi;

            /* reconstruction of each element */
            for(n=0;n<Np;n++){
                real dx = (real)(region->x[k][n] - xc);
                real dy = (real)(region->y[k][n] - yc);

                f_Q[n*Nfield+fld] = qmean + dx*qx + dy*qy;
            }

        }
    }
#if DEBUG
    fclose(fp);
#endif
    return;
}

//void SLLoc2d(PhysDomain2d *phys, double beta){
//
//    MultiReg2d   *mesh  = phys->mesh;
//    StdRegions2d *shape = mesh->stdcell;
//
//    const int Nfields   = phys->Nfields;
//    const int Np        = shape->Np;
//    const int K         = mesh->K;
//    const int Nfaces    = shape->Nfaces;
//    const int Nfp       = shape->Nfp;
//    const int nprocs    = mesh->nprocs;
//    real   *f_Q  = phys->f_Q;
//    double *area = mesh->area;
//
//    int fld,i,k,f;
//
//    /* fetch boundary values */
//    MPI_Request *mpi_out_requests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
//    MPI_Request *mpi_in_requests  = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
//    MPI_Status  *outstatus        = (MPI_Status* ) calloc(nprocs, sizeof(MPI_Status));
//    MPI_Status  *instatus         = (MPI_Status*) calloc(nprocs, sizeof(MPI_Status));
//
//    real *qmean = (real *) calloc(K, sizeof(real));
//    real *qpx   = (real *) calloc(K, sizeof(real)); /* unlimited gradient on x */
//    real *qpy   = (real *) calloc(K, sizeof(real)); /* unlimited gradient on y */
//    real *qin   = (real *) calloc(mesh->parEtotalout, sizeof(real));
//    real *qout  = (real *) calloc(mesh->parEtotalout, sizeof(real));
//
//    for (fld=0;fld<Nfields;fld++) {
//
//        /* evaluate cell mean value */
//        int sj = 0, ind;
//        for (k = 0; k < K; k++) {
//            qmean[k] = 0.0;
//            for (i=0;i<Np;i++) {
//                ind = (k*Np+i)*Nfields + fld;
//                qmean[k] += shape->wv[i] * mesh->J[sj++] * f_Q[ind];
//            }
//            qmean[k] /= area[k];
//        }
//
//        int Nmess;
//        FetchParmapEle2d(phys, qmean, qin, qout, mpi_out_requests, mpi_in_requests, &Nmess);
//
//        /* get the unlimited gradient */
//        for (k=0; k<K; k++) {
//            qpx[k] = 0.0;
//            qpy[k] = 0.0;
//
//            for (f = 0; f < Nfaces; f++) {
//                /* vertex index */
//                int l1 = shape->Fmask[f][0];       /* local element indices */
//                int l2 = shape->Fmask[f][Nfp - 1]; /* local element indices */
//                int v1 = (k * Np + l1) * Nfields + fld;
//                int v2 = (k * Np + l2) * Nfields + fld;
//                /* mean value on edge */
//                real qmean1 = (real) ((f_Q[v1] + f_Q[v2]) * 0.5);
//
//                real dx = (real) (mesh->x[k][l2] - mesh->x[k][l1]);
//                real dy = (real) (mesh->y[k][l2] - mesh->y[k][l1]);
//
//                qpx[k] += qmean1 * dy;
//                qpy[k] -= qmean1 * dx;
//            }
//            qpx[k] /= area[k];
//            qpy[k] /= area[k];
//        }
//
//        MPI_Waitall(Nmess, mpi_out_requests, outstatus);
//        MPI_Waitall(Nmess, mpi_in_requests, instatus);
//
//        /* compute the limited gradient */
//        sj = 0;
//        for(k=0;k<K;k++) {
//            real qmax = qmean[k];
//            real qmin = qmean[k];
//            real qnear; /* mean value of adjacent element */
//            for (f=0;f<Nfaces;f++) {
//                int idP = mesh->EToE[k][f];
//                if(idP>=0){
//                    qnear = qmean[idP];
//                }else{
//                    idP = (-1-idP);
//                    qnear = qin[idP];
//                }
//                qmax = max(qmax, qnear);
//                qmin = min(qmin, qnear);
//            }
//
//            real xc=0,yc=0;
//            double t, psi = 1.0;
//            for(i=0;i<Np;i++){
//                ind = (k*Np + i)*Nfields + fld;
//                real qval = f_Q[ind], rk = 1.0;
//                /* compute the limiter `psi` */
//                if(qval > qmax){
//                    rk = (qmax - qmean[k])/(qval - qmean[k]);
//                }else if(qval < qmin){
//                    rk = (qmin - qmean[k])/(qval - qmean[k]);
//                }
//                t   = max(min(beta*rk, 1.0), min(rk, beta));
//                psi = min(psi, t);
//
//                /* compute the centre of cell */
//                xc += shape->wv[i]*mesh->J[sj  ]*mesh->x[k][i];
//                yc += shape->wv[i]*mesh->J[sj++]*mesh->y[k][i];
//            }
//            xc /= area[k];
//            yc /= area[k];
//
//            /* compute the limited gradient */
//            qpx[k] *= psi;
//            qpy[k] *= psi;
//
//            /* reconstruction of each element */
//            for(i=0;i<Np;i++){
//                ind = (k*Np + i)*Nfields + fld;
//                real dx = (real)(mesh->x[k][i] - xc);
//                real dy = (real)(mesh->y[k][i] - yc);
//
//                f_Q[ind] = qmean[k] + dx*qpx[k] + dy*qpy[k];
//            }
//        }
//    }
//    /* deallocation */
//    free(qmean);
//    free(qpx);  free(qpy);
//    free(qin);  free(qout);
//
//    free(outstatus);
//    free(mpi_out_requests);
//    free(mpi_in_requests);
//}


/**
 * @brief
 * Slope limiter from Anastasiou and Chan (1997) for two dimensional problems.
 *
 * @details
 * The slope limiter will act on each physical variables and reconstruct the
 * scalar distribution. The difference between `SL2d` and `SLLoc2d` lays on the
 * calculation of boundary values. In `SL2d`, the boundary value is evaluated by
 * the inverse distance weighting method, which is given by
 * \f[ h_b = w_1 h_1 + w_2 h_2 \f]
 * where \f$ w_1 = d_2/(d_1+d_2), \, w_2 = d_1/(d_1+d_2) \f$, and
 * \f$ d_1 = \mathbf{x_1 - x_b}^2, \, d_2 = \mathbf{x_2 - x_b}^2  \f$.
 *
 * @param[in] phys
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * phys   |   PhysDomain2d*  | contains reconstruct values
 *
 */
//void SL2d(PhysDomain2d *phys, double beta){
//
//    MultiReg2d   *mesh  = phys->mesh;
//    StdRegions2d *shape = mesh->stdcell;
//
//    const int Nfields   = phys->Nfields;
//    const int Np        = shape->Np;
//    const int K         = mesh->K;
//    const int Nfaces    = shape->Nfaces;
//    const int Nfp       = shape->Nfp;
//    const int nprocs    = mesh->nprocs;
//    real   *f_Q  = phys->f_Q;
//    double *area = mesh->area;
//
//    int fld,i,k,f;
//
//    /* fetch boundary values */
//    MPI_Request *mpi_out_requests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
//    MPI_Request *mpi_in_requests  = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
//    MPI_Status  *outstatus        = (MPI_Status* ) calloc(nprocs, sizeof(MPI_Status));
//    MPI_Status  *instatus         = (MPI_Status*) calloc(nprocs, sizeof(MPI_Status));
//
//    real *qmean = (real *) calloc(K, sizeof(real));
//    real *qpx   = (real *) calloc(K, sizeof(real)); /* unlimited gradient on x */
//    real *qpy   = (real *) calloc(K, sizeof(real)); /* unlimited gradient on y */
//    real *qin   = (real *) calloc(mesh->parEtotalout, sizeof(real));
//    real *qout  = (real *) calloc(mesh->parEtotalout, sizeof(real));
//    real *xc    = (real *) calloc(K, sizeof(real));
//    real *xcin  = (real *) calloc(mesh->parEtotalout, sizeof(real));
//    real *xcout = (real *) calloc(mesh->parEtotalout, sizeof(real));
//    real *yc    = (real *) calloc(K, sizeof(real));
//    real *ycin  = (real *) calloc(mesh->parEtotalout, sizeof(real));
//    real *ycout = (real *) calloc(mesh->parEtotalout, sizeof(real));
//
//    for (fld=0;fld<Nfields;fld++) {
//
//        /* evaluate cell mean value */
//        int sj = 0, ind;
//        for (k = 0; k < K; k++) {
//            qmean[k] = 0.0;
//            for (i=0;i<Np;i++) {
//                ind = (k*Np+i)*Nfields + fld;
//                /* compute the mean value of scalar field */
//                qmean[k] += shape->wv[i] * mesh->J[sj] * f_Q[ind];
//
//                /* compute the centre of cell */
//                xc[k] += shape->wv[i]*mesh->J[sj  ]*mesh->x[k][i];
//                yc[k] += shape->wv[i]*mesh->J[sj++]*mesh->y[k][i];
//            }
//            qmean[k] /= area[k];
//            xc[k] /= area[k];
//            yc[k] /= area[k];
//        }
//
//        int Nmess;
//        FetchParmapEle2d(phys, xc, xcin, xcout, mpi_out_requests, mpi_in_requests, &Nmess);
//        MPI_Waitall(Nmess, mpi_out_requests, outstatus);
//        MPI_Waitall(Nmess, mpi_in_requests, instatus);
//
//        FetchParmapEle2d(phys, yc, ycin, ycout, mpi_out_requests, mpi_in_requests, &Nmess);
//        MPI_Waitall(Nmess, mpi_out_requests, outstatus);
//        MPI_Waitall(Nmess, mpi_in_requests, instatus);
//
//        FetchParmapEle2d(phys, qmean, qin, qout, mpi_out_requests, mpi_in_requests, &Nmess);
//        MPI_Waitall(Nmess, mpi_out_requests, outstatus);
//        MPI_Waitall(Nmess, mpi_in_requests, instatus);
//
//        /* get the unlimited gradient */
//        for (k=0; k<K; k++) {
//            qpx[k] = 0.0;
//            qpy[k] = 0.0;
//
//            for (f = 0; f < Nfaces; f++) {
//                int td = mesh->EToE[k][f];  /* index of adjacent element */
//                int l1 = shape->Fmask[f][0];
//                int l2 = shape->Fmask[f][Nfp - 1];
//
//                /* mean value on edge */
//                real d1, d2, q1, q2;
//                real xb, yb, x1, y1, x2, y2;
//                xb = (real) ((mesh->x[k][l2] + mesh->x[k][l1])*0.5);
//                yb = (real) ((mesh->y[k][l2] + mesh->y[k][l1])*0.5);
//                x1 = (real) xc[k];
//                y1 = (real) yc[k];
//                q1 = qmean[k];
//                if(td>=0){
//                    x2 = xc[td];
//                    y2 = yc[td];
//                    q2 = qmean[td];
//                }else{
//                    td = (-1-td);
//                    x2 = xcin[td];
//                    y2 = ycin[td];
//                    q2 = qin [td];
//                }
//
//                d1 = (x1 - xb)*(x1 - xb) + (y1 - yb)*(y1 - yb);
//                d2 = (x2 - xb)*(x2 - xb) + (y2 - yb)*(y2 - yb);
//
//                real qmean1 = (q1*d2 + q2*d1)/(d1 + d2);
//
//                real dx = (real) (mesh->x[k][l2] - mesh->x[k][l1]);
//                real dy = (real) (mesh->y[k][l2] - mesh->y[k][l1]);
//
//                qpx[k] += qmean1 * dy;
//                qpy[k] -= qmean1 * dx;
//            }
//            qpx[k] /= area[k];
//            qpy[k] /= area[k];
//        }
//
//        /* compute the limited gradient */
//        for(k=0;k<K;k++) {
//            real qmax = qmean[k];
//            real qmin = qmean[k];
//            real qnear; /* mean value of adjacent element */
//            for (f=0;f<Nfaces;f++) {
//                int idP = mesh->EToE[k][f];
//                if(idP>=0){
//                    qnear = qmean[idP];
//                }else{
//                    idP = (-1-idP);
//                    qnear = qin[idP];
//                }
//                qmax = max(qmax, qnear);
//                qmin = min(qmin, qnear);
//            }
//
//            double t, psi = 1.0;
//            for(i=0;i<Np;i++){
//                ind = (k*Np + i)*Nfields + fld;
//                real qval = f_Q[ind], rk = 1.0;
//                /* compute the limiter `psi` */
//                if(qval > qmax){
//                    rk = (qmax - qmean[k])/(qval - qmean[k]);
//                }else if(qval < qmin){
//                    rk = (qmin - qmean[k])/(qval - qmean[k]);
//                }
//                t   = max(min(beta*rk, 1.0), min(rk, beta));
//                psi = min(psi, t);
//            }
//
//            /* compute the limited gradient */
//            qpx[k] *= psi;
//            qpy[k] *= psi;
//
//            /* reconstruction */
//            for(i=0;i<Np;i++){
//                ind = (k*Np + i)*Nfields + fld;
//                real dx = (real)(mesh->x[k][i] - xc[k]);
//                real dy = (real)(mesh->y[k][i] - yc[k]);
//
//                f_Q[ind] = qmean[k] + dx*qpx[k] + dy*qpy[k];
//            }
//        }
//    }
//
//    free(qmean);
//    free(qpx);  free(qpy);
//    free(qin);  free(qout);
//
//    free(outstatus);
//    free(mpi_out_requests);
//    free(mpi_in_requests);
//}