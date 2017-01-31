//
// Created by li12242 on 16/8/3.
//

#include <MultiRegions/mr_mesh.h>
#include "pf_limiter.h"
#include "pf_cellMean.h"
#include "pf_fetchBuffer.h"

#define DEBUG 0
#if DEBUG
#include "LibUtilities/UTest.h"
#endif

#define TOL 1e-8
#define EXCEPTION(u, u1, u2) ( ((u-u1)>TOL)& ((u-u2)>TOL) )||( ((u1-u)>TOL)&((u2-u)>TOL))


static void pf_faceMean(physField *phys, int k, int f, real *f_mean){

    const int Nfield = phys->Nfield;
    const int Nfp = phys->cell->Nfp;
    const int Np = phys->cell->Np;

    double *ws = phys->cell->ws;
    real *f_Q = phys->f_Q + k*Np*Nfield;
    int *fmask = phys->cell->Fmask[f];

    register int n,m,fld;
    /* initialization */
    for(fld=0;fld<Nfield;fld++){
        f_mean[fld] = 0;
    }
    for(n=0;n<Nfp;n++){
        m = fmask[n]*Nfield;
        double w = ws[n]*0.5;
        for(fld=0;fld<Nfield;fld++){
            f_mean[fld] += f_Q[m+fld]*w;
        }
    }
}

/**
 * @brief calculate the gradient by Gauss-Green method
 *
 * @param[in] phys
 * @param[in] uf_mean
 * @param[out] px
 * @param[out] py
 */
static void phys_gradient_Green(physField *phys, real *px, real *py){
    const int K = phys->grid->K;
    const int Nfield = phys->Nfield;
    const int Nfaces = phys->cell->Nfaces;
    const int Nfp = phys->cell->Nfp;

    double **x = phys->region->x;
    double **y = phys->region->y;

    register int k,f,fld,sk;
    for(k=0;k<K;k++){
        double A = 1.0/phys->region->size[k];

        for(fld=0;fld<Nfield;fld++){
            sk = k*Nfield+fld;
            px[sk] = 0;
            py[sk] = 0;
        }

        for(f=0;f<Nfaces;f++){
            real f_mean[Nfield];
            pf_faceMean(phys, k, f, f_mean);

            int v1 = phys->cell->Fmask[f][0];
            int v2 = phys->cell->Fmask[f][Nfp-1];
            double x1 = x[k][v1];
            double x2 = x[k][v2];
            double y1 = y[k][v1];
            double y2 = y[k][v2];

            double dx =  (y2-y1);
            double dy = -(x2-x1);

            for(fld=0;fld<Nfield;fld++){
                sk = k*Nfield+fld;
                px[sk] += f_mean[fld]*dx*A;
                py[sk] += f_mean[fld]*dy*A;
#if DEBUG
                if(!phys->mesh->procid)
                    printf("k=%d, f=%d, fld=%d, f_mean=%f, px=%f, py=%f\n",
                           k,f,fld,f_mean[fld],px[sk],py[sk]);
#endif
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
void pf_slloc2d(physField *phys, double beta){

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
    pf_cellMean(phys);

    /* 2. fetch cell info with other processes */
    MPI_Request mpi_out_requests[nprocs];
    MPI_Request mpi_in_requests[nprocs];
    int Nmess;

    /* do sends and recv */
    pf_fetchCellBuffer(phys, mpi_out_requests, mpi_in_requests, &Nmess);

    /* 4.1. judge the trouble cell, get c_max and c_min  */
    int **EToE = phys->mesh->EToE;
    int **EToP = phys->mesh->EToP;
    real cell_max[K*Nfield], cell_min[K*Nfield];
    int t_cell[K];

    for(k=0;k<K;k++){
        t_cell[k] = 0;
        real *c_mean = phys->c_Q + k*Nfield; // mean value of k-th cell
        real c_max[Nfield];
        real c_min[Nfield];
        /* initialize of c_max & c_min */
        for(fld=0;fld<Nfield;fld++){
            c_max[fld] = c_mean[fld];
            c_min[fld] = c_mean[fld];
        }

        for(f=0;f<Nfaces;f++){
            int e = EToE[k][f];
            int p = EToP[k][f];
            // if adjacent cell is on other process, jump this cycle
            if(p!=procid) continue;

            real f_mean[Nfield];
            pf_faceMean(phys, k, f, f_mean);

//            if(e!=k){
            for(fld=0;fld<Nfield;fld++){
                real c_next = phys->c_Q[e*Nfield+fld];

                c_max[fld] = max(c_max[fld], c_next);
                c_min[fld] = min(c_min[fld], c_next);

                if( EXCEPTION(f_mean[fld], c_mean[fld], c_next) )
                    t_cell[k] = 1;
#if DEBUG
                    fprintf(fp, "k=%d, f=%d, fld=%d, c_mean=%f, c_next=%f, f_mean=%f\n",
                            k, f, fld, c_mean[fld], c_next, f_mean[fld]);
#endif
            }
//            }else if(e==k){
//                for(fld=0;fld<Nfield;fld++){
//                    c_max[fld] = max(c_max[fld], f_mean[fld]);
//                    c_min[fld] = min(c_min[fld], f_mean[fld]);
//                }
//            }
        }

        int sk = k*Nfield;
        for(fld=0;fld<Nfield;fld++){
            cell_max[sk+fld] = c_max[fld];
            cell_min[sk+fld] = c_min[fld];
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
        real f_mean[Nfield];
        pf_faceMean(phys, k, f, f_mean);
        for(fld=0;fld<Nfield;fld++){
            real c_next = phys->c_inQ[n*Nfield+fld];
            int sk = k*Nfield+fld;
            cell_max[sk] = max(cell_max[sk], c_next);
            cell_min[sk] = min(cell_min[sk], c_next);
            if( EXCEPTION(f_mean[fld], c_mean[fld], c_next) ){
#if DEBUG
                fprintf(fp, "k=%d, f=%d, fld=%d, c_mean=%f, c_next=%f, f_mean=%f\n",
                        k, f, fld, c_mean[fld], c_next, f_mean[fld]);
#endif
                t_cell[k] = 1;
            }
        }
    }

#if DEBUG
    PrintIntVector2File(fp, "trouble_cell", t_cell, K);
#endif

    /* 5. calculate the unlimited gradient */
    real px[K*Nfield], py[K*Nfield];
    phys_gradient_Green(phys, px, py);
#if DEBUG
    PrintVector2File(fp, "px", px, K*Nfield);
    PrintVector2File(fp, "py", py, K*Nfield);
#endif

    /* 6. calculate the limited results */
    multiReg *region = phys->region;
    for(k=0;k<K;k++){
        /* check if the cell need limited */
        //if(t_cell[k]==0 ) continue;

        double A = 1.0/phys->region->size[k];

        double xc = A * mr_reg_integral(region, k, region->x[k]);
        double yc = A * mr_reg_integral(region, k, region->y[k]);

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
#if DEBUG
                fprintf(fp, "k=%d, fld=%d, n=%d, rk=%f, psi=%f\n",k,fld,n,rk,psi);
#endif
            }

            /* compute the limited gradient */
            real qx = px[k*Nfield+fld] * psi;
            real qy = py[k*Nfield+fld] * psi;

            /* reconstruction of each element */
            for(n=0;n<Np;n++){
                real dx = (real)(region->x[k][n] - xc);
                real dy = (real)(region->y[k][n] - yc);

                f_Q[n*Nfield+fld] = qmean + dx*qx + dy*qy;
#if DEBUG
                fprintf(fp, "k=%d, fld=%d, n=%d, px=%f, py=%f, dx=%f, dy=%f, f=%f\n",
                        k,fld,n,qx,qy,dx,dy,f_Q[n*Nfield+fld]);
#endif
            }

        }
    }
#if 0
    k=1108; n=1, fld=1;
    printf("procid=%d, qmax=%f, qmin=%f, qmean=%f\n", phys->mesh->procid,
           cell_max[k*Nfield+fld], cell_min[k*Nfield+fld], phys->c_Q[k*Nfield+fld]);
    int ind = 13302-3;
    printf("procid=%d, var=[%f, %f, %f], var=[%f, %f, %f]\n",
           phys->mesh->procid,
           phys->f_Q[ind  ], phys->f_Q[ind+1], phys->f_Q[ind+2],
           phys->f_Q[ind+3], phys->f_Q[ind+4], phys->f_Q[ind+5]);
#endif
#if DEBUG
    fclose(fp);
#endif
    return;
}