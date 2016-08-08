//
// Created by li12242 on 16/8/3.
//

#include "SlopeLimiter.h"

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
void SLLoc2d(PhysDomain2d *phys, double beta){

    MultiReg2d   *mesh  = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    const int Nfields   = phys->Nfields;
    const int Np        = shape->Np;
    const int K         = mesh->K;
    const int Nfaces    = shape->Nfaces;
    const int Nfp       = shape->Nfp;
    const int nprocs    = mesh->nprocs;
    real   *f_Q  = phys->f_Q;
    double *area = mesh->area;

    int fld,i,k,f;

    /* fetch boundary values */
    MPI_Request *mpi_out_requests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_in_requests  = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Status  *outstatus        = (MPI_Status* ) calloc(nprocs, sizeof(MPI_Status));
    MPI_Status  *instatus         = (MPI_Status*) calloc(nprocs, sizeof(MPI_Status));

    real *qmean = (real *) calloc(K, sizeof(real));
    real *qpx   = (real *) calloc(K, sizeof(real)); /* unlimited gradient on x */
    real *qpy   = (real *) calloc(K, sizeof(real)); /* unlimited gradient on y */
    real *qin   = (real *) calloc(mesh->parEtotalout, sizeof(real));
    real *qout  = (real *) calloc(mesh->parEtotalout, sizeof(real));

    for (fld=0;fld<Nfields;fld++) {

        /* evaluate cell mean value */
        int sj = 0, ind;
        for (k = 0; k < K; k++) {
            qmean[k] = 0.0;
            for (i=0;i<Np;i++) {
                ind = (k*Np+i)*Nfields + fld;
                qmean[k] += shape->wv[i] * mesh->J[sj++] * f_Q[ind];
            }
            qmean[k] /= area[k];
        }

        /* get the unlimited gradient */
        for (k=0; k<K; k++) {
            qpx[k] = 0.0;
            qpy[k] = 0.0;

            for (f = 0; f < Nfaces; f++) {
                /* vertex index */
                int l1 = shape->Fmask[f][0];
                int l2 = shape->Fmask[f][Nfp - 1];
                int v1 = (k * Np + l1) * Nfields + fld;
                int v2 = (k * Np + l2) * Nfields + fld;
                /* mean value on edge */
                real qmean1 = (real) ((f_Q[v1] + f_Q[v2]) * 0.5);

                real dx = (real) (mesh->x[k][l2] - mesh->x[k][l1]);
                real dy = (real) (mesh->y[k][l2] - mesh->y[k][l1]);

                qpx[k] += qmean1 * dy;
                qpy[k] -= qmean1 * dx;
            }
            qpx[k] /= area[k];
            qpy[k] /= area[k];
        }

        int Nmess;
        FetchParmapEle2d(phys, qmean, qin, qout, mpi_out_requests, mpi_in_requests, &Nmess);
        MPI_Waitall(Nmess, mpi_out_requests, outstatus);
        MPI_Waitall(Nmess, mpi_in_requests, instatus);

        /* compute the limited gradient */
        sj = 0;
        for(k=0;k<K;k++) {
            real qmax = qmean[k];
            real qmin = qmean[k];
            real qnear; /* mean value of adjacent element */
            for (f=0;f<Nfaces;f++) {
                int idP = mesh->EToE[k][f];
                if(idP>=0){
                    qnear = qmean[idP];
                }else{
                    idP = (-1-idP);
                    qnear = qin[idP];
                }
                qmax = max(qmax, qnear);
                qmin = min(qmin, qnear);
            }

            real xc=0,yc=0;
            double t, psi = 1.0;
            for(i=0;i<Np;i++){
                ind = (k*Np + i)*Nfields + fld;
                real qval = f_Q[ind], rk = 1.0;
                /* compute the limiter `psi` */
                if(qval > qmax){
                    rk = (qmax - qmean[k])/(qval - qmean[k]);
                }else if(qval < qmin){
                    rk = (qmin - qmean[k])/(qval - qmean[k]);
                }
                t   = max(min(beta*rk, 1.0), min(rk, beta));
                psi = min(psi, t);

                /* compute the centre of cell */
                xc += shape->wv[i]*mesh->J[sj  ]*mesh->x[k][i];
                yc += shape->wv[i]*mesh->J[sj++]*mesh->y[k][i];
            }
            xc /= area[k];
            yc /= area[k];

            /* compute the limited gradient */
            qpx[k] *= psi;
            qpy[k] *= psi;

            /* reconstruction of each element */
            for(i=0;i<Np;i++){
                ind = (k*Np + i)*Nfields + fld;
                real dx = (real)(mesh->x[k][i] - xc);
                real dy = (real)(mesh->y[k][i] - yc);

                f_Q[ind] = qmean[k] + dx*qpx[k] + dy*qpy[k];
            }
        }
    }
    /* deallocation */
    free(qmean);
    free(qpx);  free(qpy);
    free(qin);  free(qout);

    free(outstatus);
    free(mpi_out_requests);
    free(mpi_in_requests);
}


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
void SL2d(PhysDomain2d *phys, double beta){

    MultiReg2d   *mesh  = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    const int Nfields   = phys->Nfields;
    const int Np        = shape->Np;
    const int K         = mesh->K;
    const int Nfaces    = shape->Nfaces;
    const int Nfp       = shape->Nfp;
    const int nprocs    = mesh->nprocs;
    real   *f_Q  = phys->f_Q;
    double *area = mesh->area;

    int fld,i,k,f;

    /* fetch boundary values */
    MPI_Request *mpi_out_requests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_in_requests  = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Status  *outstatus        = (MPI_Status* ) calloc(nprocs, sizeof(MPI_Status));
    MPI_Status  *instatus         = (MPI_Status*) calloc(nprocs, sizeof(MPI_Status));

    real *qmean = (real *) calloc(K, sizeof(real));
    real *qpx   = (real *) calloc(K, sizeof(real)); /* unlimited gradient on x */
    real *qpy   = (real *) calloc(K, sizeof(real)); /* unlimited gradient on y */
    real *qin   = (real *) calloc(mesh->parEtotalout, sizeof(real));
    real *qout  = (real *) calloc(mesh->parEtotalout, sizeof(real));
    real *xc    = (real *) calloc(K, sizeof(real));
    real *xcin  = (real *) calloc(mesh->parEtotalout, sizeof(real));
    real *xcout = (real *) calloc(mesh->parEtotalout, sizeof(real));
    real *yc    = (real *) calloc(K, sizeof(real));
    real *ycin  = (real *) calloc(mesh->parEtotalout, sizeof(real));
    real *ycout = (real *) calloc(mesh->parEtotalout, sizeof(real));

    for (fld=0;fld<Nfields;fld++) {

        /* evaluate cell mean value */
        int sj = 0, ind;
        for (k = 0; k < K; k++) {
            qmean[k] = 0.0;
            for (i=0;i<Np;i++) {
                ind = (k*Np+i)*Nfields + fld;
                /* compute the mean value of scalar field */
                qmean[k] += shape->wv[i] * mesh->J[sj] * f_Q[ind];

                /* compute the centre of cell */
                xc[k] += shape->wv[i]*mesh->J[sj  ]*mesh->x[k][i];
                yc[k] += shape->wv[i]*mesh->J[sj++]*mesh->y[k][i];
            }
            qmean[k] /= area[k];
            xc[k] /= area[k];
            yc[k] /= area[k];
        }

        int Nmess;
        FetchParmapEle2d(phys, xc, xcin, xcout, mpi_out_requests, mpi_in_requests, &Nmess);
        MPI_Waitall(Nmess, mpi_out_requests, outstatus);
        MPI_Waitall(Nmess, mpi_in_requests, instatus);

        FetchParmapEle2d(phys, yc, ycin, ycout, mpi_out_requests, mpi_in_requests, &Nmess);
        MPI_Waitall(Nmess, mpi_out_requests, outstatus);
        MPI_Waitall(Nmess, mpi_in_requests, instatus);

        FetchParmapEle2d(phys, qmean, qin, qout, mpi_out_requests, mpi_in_requests, &Nmess);
        MPI_Waitall(Nmess, mpi_out_requests, outstatus);
        MPI_Waitall(Nmess, mpi_in_requests, instatus);

        /* get the unlimited gradient */
        for (k=0; k<K; k++) {
            qpx[k] = 0.0;
            qpy[k] = 0.0;

            for (f = 0; f < Nfaces; f++) {
                int td = mesh->EToE[k][f];  /* index of adjacent element */
                int l1 = shape->Fmask[f][0];
                int l2 = shape->Fmask[f][Nfp - 1];

                /* mean value on edge */
                real d1, d2, q1, q2;
                real xb, yb, x1, y1, x2, y2;
                xb = (real) ((mesh->x[k][l2] + mesh->x[k][l1])*0.5);
                yb = (real) ((mesh->y[k][l2] + mesh->y[k][l1])*0.5);
                x1 = (real) xc[k];
                y1 = (real) yc[k];
                q1 = qmean[k];
                if(td>=0){
                    x2 = xc[td];
                    y2 = yc[td];
                    q2 = qmean[td];
                }else{
                    td = (-1-td);
                    x2 = xcin[td];
                    y2 = ycin[td];
                    q2 = qin [td];
                }

                d1 = (x1 - xb)*(x1 - xb) + (y1 - yb)*(y1 - yb);
                d2 = (x2 - xb)*(x2 - xb) + (y2 - yb)*(y2 - yb);

                real qmean1 = (q1*d2 + q2*d1)/(d1 + d2);

                real dx = (real) (mesh->x[k][l2] - mesh->x[k][l1]);
                real dy = (real) (mesh->y[k][l2] - mesh->y[k][l1]);

                qpx[k] += qmean1 * dy;
                qpy[k] -= qmean1 * dx;
            }
            qpx[k] /= area[k];
            qpy[k] /= area[k];
        }

        /* compute the limited gradient */
        for(k=0;k<K;k++) {
            real qmax = qmean[k];
            real qmin = qmean[k];
            real qnear; /* mean value of adjacent element */
            for (f=0;f<Nfaces;f++) {
                int idP = mesh->EToE[k][f];
                if(idP>=0){
                    qnear = qmean[idP];
                }else{
                    idP = (-1-idP);
                    qnear = qin[idP];
                }
                qmax = max(qmax, qnear);
                qmin = min(qmin, qnear);
            }

            double t, psi = 1.0;
            for(i=0;i<Np;i++){
                ind = (k*Np + i)*Nfields + fld;
                real qval = f_Q[ind], rk = 1.0;
                /* compute the limiter `psi` */
                if(qval > qmax){
                    rk = (qmax - qmean[k])/(qval - qmean[k]);
                }else if(qval < qmin){
                    rk = (qmin - qmean[k])/(qval - qmean[k]);
                }
                t   = max(min(beta*rk, 1.0), min(rk, beta));
                psi = min(psi, t);
            }

            /* compute the limited gradient */
            qpx[k] *= psi;
            qpy[k] *= psi;

            /* reconstruction */
            for(i=0;i<Np;i++){
                ind = (k*Np + i)*Nfields + fld;
                real dx = (real)(mesh->x[k][i] - xc[k]);
                real dy = (real)(mesh->y[k][i] - yc[k]);

                f_Q[ind] = qmean[k] + dx*qpx[k] + dy*qpy[k];
            }
        }
    }

    free(qmean);
    free(qpx);  free(qpy);
    free(qin);  free(qout);

    free(outstatus);
    free(mpi_out_requests);
    free(mpi_in_requests);
}