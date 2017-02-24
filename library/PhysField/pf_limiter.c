//
// Created by li12242 on 16/8/3.
//

#include <MultiRegions/mr_mesh.h>
#include "pf_limiter.h"
#include "pf_cellMean.h"
#include "pf_fetchBuffer.h"

#define DEBUG 0
#if DEBUG
#include "Utility/UTest.h"
#endif

#define TOL 1e-8
#define INFTY 1e10
#define EXCEPTION(u, u1, u2) ( ((u-u1)>TOL)& ((u-u2)>TOL) )||( ((u1-u)>TOL)&((u2-u)>TOL))

/**
 * @brief integral averaged values of each faces
 * @param[in] phys physical field structure.
 * @param[in,out] f_mean integral averaged values of each faces
 */
static void pf_face_mean(physField *phys, dg_real *f_mean){

    const int Nfield = phys->Nfield;
    const int Nfaces = phys->cell->Nfaces;
    const int Nfp = phys->cell->Nfp;
    const int Np = phys->cell->Np;
    const int K = phys->grid->K;

    double *ws = phys->cell->ws;
    register int k,f,n,m,fld;
    for(k=0;k<K;k++){
        dg_real *f_Q = phys->f_Q + k*Np*Nfield;
        for(f=0;f<Nfaces;f++){
            int sk = k*Nfaces*Nfield + f*Nfield;
            /* initialization */
            for(fld=0;fld<Nfield;fld++){
                f_mean[sk+fld] = 0;
            }
            int *fmask = phys->cell->Fmask[f];
            for(n=0;n<Nfp;n++){
                m = fmask[n]*Nfield;
                double w = ws[n]*0.5;
                for(fld=0;fld<Nfield;fld++){
                    f_mean[sk+fld] += f_Q[m+fld]*w;
                }
            }
        }
    }
    return;
}

static void pf_weiface_mean(physField *phys, dg_real *f_mean){
    const int Nfield = phys->Nfield;
    const int Nfaces = phys->cell->Nfaces;
    const int K = phys->grid->K;
    const int Nfp = phys->cell->Nfp;
    const int procid = phys->mesh->procid;
    const int nprocs = phys->mesh->nprocs;

    multiReg *region = phys->region;
    parallMesh *mesh = phys->mesh;
    register int k,f,n,fld;
    dg_real *xc = vector_real_create(K);
    dg_real *yc = vector_real_create(K);
    for(k=0;k<K;k++){
        const double Area = 1.0/region->size[k];
        xc[k] = Area * mr_reg_integral(region, k, region->x[k]);
        yc[k] = Area * mr_reg_integral(region, k, region->y[k]);
    }
    const int parallCellNum = phys->mesh->parallCellNum;
    dg_real *xc_out = vector_real_create(parallCellNum);
    dg_real *xc_in = vector_real_create(parallCellNum);
    dg_real *yc_out = vector_real_create(parallCellNum);
    dg_real *yc_in = vector_real_create(parallCellNum);

    for(n=0;n<parallCellNum;++n){
        int sk = mesh->cellIndexOut[n];
        xc_out[n] = xc[sk];
        yc_out[n] = yc[sk];
    }
    MPI_Request xc_out_requests[nprocs];
    MPI_Request xc_in_requests[nprocs];
    MPI_Request yc_out_requests[nprocs];
    MPI_Request yc_in_requests[nprocs];
    int Nmess;
    pf_fetchBuffer(procid, nprocs, mesh->Npar, xc_out, xc_in,
                   xc_out_requests, xc_in_requests, &Nmess);
    pf_fetchBuffer(procid, nprocs, mesh->Npar, yc_out, yc_in,
                   yc_out_requests, yc_in_requests, &Nmess);

    dg_real *c_Q = phys->c_Q;
    int **EToE = mesh->EToE;
    int **EToP = mesh->EToP;
    double **x = phys->region->x;
    double **y = phys->region->y;
    for(k=0;k<K;k++){
        for(f=0;f<Nfaces;f++){
            int sk = k*Nfaces*Nfield + f*Nfield;
            /* initialization */
            for(fld=0;fld<Nfield;fld++){
                f_mean[sk+fld] = 0;
            }

            int e = EToE[k][f];
            int p = EToP[k][f];
            if(p!=procid) { continue; }

            int *fmask = phys->cell->Fmask[f];
            int v1 = fmask[0];
            int v2 = fmask[Nfp-1];
            double xf = (x[k][v1] + x[k][v2])*0.5;
            double yf = (y[k][v1] + y[k][v2])*0.5;

            double d1 = (xf - xc[k])*(xf - xc[k]) + (yf - yc[k])*(yf - yc[k]);
            double d2 = (xf - xc[e])*(xf - xc[e]) + (yf - yc[e])*(yf - yc[e]);
            dg_real w1 = d2/(d1 + d2);
            dg_real w2 = d1/(d1 + d2);
            for(fld=0;fld<Nfield;fld++)
                f_mean[sk+fld] = (c_Q[k*Nfield+fld]*w1 + c_Q[e*Nfield+fld]*w2);
        }
    }
    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, xc_in_requests, instatus);
    MPI_Waitall(Nmess, xc_out_requests, instatus);
    MPI_Waitall(Nmess, yc_in_requests, instatus);
    MPI_Waitall(Nmess, yc_out_requests, instatus);

    for(n=0;n<parallCellNum;n++){
        k = mesh->cellIndexIn[n];
        f = mesh->faceIndexIn[n];

        dg_real xc_next = xc_in[n];
        dg_real yc_next = yc_in[n];
        int *fmask = phys->cell->Fmask[f];
        int v1 = fmask[0];
        int v2 = fmask[Nfp-1];
        double xf = (x[k][v1] + x[k][v2])*0.5;
        double yf = (y[k][v1] + y[k][v2])*0.5;

        double d1 = (xf - xc[k])*(xf - xc[k]) + (yf - yc[k])*(yf - yc[k]);
        double d2 = (xf - xc_next)*(xf - xc_next) + (yf - yc_next)*(yf - yc_next);
        dg_real w1 = d2/(d1 + d2);
        dg_real w2 = d1/(d1 + d2);

        int sk = k*Nfaces*Nfield + f*Nfield;
        for(fld=0;fld<Nfield;fld++){
            f_mean[sk+fld] = (c_Q[k*Nfield+fld]*w1 + phys->c_inQ[n*Nfield+fld]*w2);
        }
    }

    vector_real_free(xc_out);
    vector_real_free(xc_in);
    vector_real_free(yc_out);
    vector_real_free(yc_in);

    vector_real_free(xc);
    vector_real_free(yc);
    return;
}

/**
 * @brief calculate the gradient by Gauss-Green method
 * @param[in] phys
 * @param[out] px
 * @param[out] py
 */
static void phys_gradient_Green(physField *phys, dg_real *px, dg_real *py){
    const int K = phys->grid->K;
    const int Nfield = phys->Nfield;
    const int Nfaces = phys->cell->Nfaces;
    const int Nfp = phys->cell->Nfp;

    register int k,f,fld,sk;
    int **fmask = phys->cell->Fmask;

    dg_real f_mean[K*Nfaces*Nfield];
    pf_face_mean(phys, f_mean);

    for(k=0;k<K;k++){
        double A = 1.0/phys->region->size[k];
        double *x = phys->region->x[k];
        double *y = phys->region->y[k];

        for(fld=0;fld<Nfield;fld++){
            sk = k*Nfield+fld;
            px[sk] = 0;
            py[sk] = 0;
        }

        for(f=0;f<Nfaces;f++){
            int v1 = fmask[f][0];
            int v2 = fmask[f][Nfp-1];
            double x1 = x[v1];
            double x2 = x[v2];
            double y1 = y[v1];
            double y2 = y[v2];

            double dx =  (y2-y1);
            double dy = -(x2-x1);

            for(fld=0;fld<Nfield;fld++){
                sk = k*Nfield+fld;
                int sf = k*Nfaces*Nfield + f*Nfield;
                px[sk] += f_mean[sf+fld]*dx*A;
                py[sk] += f_mean[sf+fld]*dy*A;
            }
        }
    }

    return;
}


void pf_adjacent_cellinfo(physField *phys, dg_real *cell_max, dg_real *cell_min){
    const int K = phys->grid->K;
    const int Nfield = phys->Nfield;
    const int procid = phys->mesh->procid;
    const int nprocs = phys->mesh->nprocs;
    const int Nfaces = phys->cell->Nfaces;

    register int k,n,f,fld;

    MPI_Request mpi_out_requests[nprocs];
    MPI_Request mpi_in_requests[nprocs];
    int Nmess;

    /* do sends and recv */
    pf_fetchCellBuffer(phys, mpi_out_requests, mpi_in_requests, &Nmess);

    /* local cell loop */
    int **EToE = phys->mesh->EToE;
    int **EToP = phys->mesh->EToP;
    for(k=0;k<K;k++){
        dg_real *c_mean = phys->c_Q + k*Nfield; // mean value of k-th cell
        dg_real c_max[Nfield];
        dg_real c_min[Nfield];
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

            for(fld=0;fld<Nfield;fld++){
                dg_real c_next = phys->c_Q[e*Nfield+fld];
                c_max[fld] = max(c_max[fld], c_next);
                c_min[fld] = min(c_min[fld], c_next);
            }
        }
        for(fld=0;fld<Nfield;fld++){
            int sk = k*Nfield + fld;
            cell_max[sk] = c_max[fld];
            cell_min[sk] = c_min[fld];
        }
    }

    /* wait for 2. sends and recv buffers */
    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_in_requests, instatus);
    MPI_Waitall(Nmess, mpi_out_requests, instatus);

    /* parallel cell loop */
    parallMesh *mesh = phys->mesh;
    for(n=0;n<mesh->parallCellNum;n++){
        k = mesh->cellIndexIn[n];
        for(fld=0;fld<Nfield;fld++){
            dg_real c_next = phys->c_inQ[n*Nfield+fld];
            int sk = k*Nfield+fld;
            cell_max[sk] = max(cell_max[sk], c_next);
            cell_min[sk] = min(cell_min[sk], c_next);
        }
    }
}

/**
 * @brief obtain Barth-Jeperson slope limiter.
 * @param phys
 * @param cell_max
 * @param cell_min
 * @param beta
 * @param psi
 */
void pf_BJ_limiter(physField *phys, dg_real *cell_max, dg_real *cell_min,
                   double beta, dg_real *psi){

    const int K = phys->grid->K;
    const int Nfield = phys->Nfield;
    const int Np = phys->cell->Np;

    register int k,fld,n;
    for(k=0;k<K;k++){
        dg_real *f_Q = phys->f_Q + k*Np*Nfield; // variable of k-th cell
        for(fld=0;fld<Nfield;fld++){
            int sk = k*Nfield+fld;
            dg_real qmax = cell_max[sk];
            dg_real qmin = cell_min[sk];
            dg_real qmean = phys->c_Q[sk];
            psi[sk] = 1.0;
            for(n=0;n<Np;n++){
                dg_real qval = f_Q[n*Nfield+fld];
                double rk = 1.0;
                /* compute the limiter `psi` */
                if(qval > qmax){
                    rk = (qmax - qmean)/(qval - qmean);
                }else if(qval < qmin){
                    rk = (qmin - qmean)/(qval - qmean);
                }
                dg_real t = max(min(beta*rk, 1.0), min(rk, beta));
                psi[sk] = min(psi[sk], t);
            }
        }
    }
    return;
}

/**
 * @brief
 * Slope limiter from Anastasiou and Chan (1997) for two dimensional problems.
 * @details
 * The slope limiter will act on each physical variables and reconstruct the
 * scalar distribution.
 * @param[in] phys
 */
void pf_slloc2d(physField *phys, double beta){

    const int K = phys->grid->K;
    const int Nfield = phys->Nfield;
    const int Np = phys->cell->Np;
    register int k,n,fld;

#if DEBUG
    FILE *fp = CreateLog("phys_slloc2d", phys->mesh->procid, phys->mesh->nprocs);
#endif

    /* 1. calculate the cell average value */
    pf_cellMean(phys);

    /* 2. fetch cell info with other processes */
    dg_real cell_max[K*Nfield], cell_min[K*Nfield];
    pf_adjacent_cellinfo(phys, cell_max, cell_min);

    /* 3. calculate the unlimited gradient */
    dg_real px[K*Nfield], py[K*Nfield];
    phys_gradient_Green(phys, px, py);

    /* 4. calculate the limited results */
    dg_real psi[K*Nfield];
    pf_BJ_limiter(phys, cell_max, cell_min, beta, psi);

    /* 5. reconstruction */
    multiReg *region = phys->region;
    for(k=0;k<K;k++){
        double A = 1.0/phys->region->size[k];
        double xc = A * mr_reg_integral(region, k, region->x[k]);
        double yc = A * mr_reg_integral(region, k, region->y[k]);

        dg_real *f_Q = phys->f_Q + k*Np*Nfield; // variable of k-th cell
        for(fld=0;fld<Nfield;fld++){
            int sk = k*Nfield+fld;
            /* compute the limited gradient */
            dg_real qx = px[sk] * psi[sk];
            dg_real qy = py[sk] * psi[sk];

            /* reconstruction of each element */
            dg_real qmean = phys->c_Q[sk];
            for(n=0;n<Np;n++){
                dg_real dx = (dg_real)(region->x[k][n] - xc);
                dg_real dy = (dg_real)(region->y[k][n] - yc);

                f_Q[n*Nfield+fld] = qmean + dx*qx + dy*qy;
            }
        }
    }
#if DEBUG
    fclose(fp);
#endif
    return;
}


static void triGradient(dg_real *x, dg_real *y, dg_real *f, dg_real *px, dg_real *py){
    dg_real a[4], b[2];
    a[0] = x[1] - x[0]; a[1] = y[1] - y[0];
    a[2] = x[2] - x[0]; a[3] = y[2] - y[0];
    b[0] = f[1] - f[0]; b[1] = f[2] - f[0];

    dg_real det = a[0]*a[3] - a[1]*a[2];
    *px = ( b[0]*a[3] - b[1]*a[1])/det;
    *py = (-b[0]*a[2] + b[1]*a[0])/det;
}

static void VA_weigrad(
        const int Nfaces, const int Nfield,
        dg_real *px, dg_real *py,
        dg_real *wpx, dg_real *wpy){

    register int i,j,fld;
    dg_real EPSILON = 1.0e-12;
    dg_real gra_x[Nfaces], gra_y[Nfaces], gra_det[Nfaces];
    for(fld=0;fld<Nfield;fld++){
        dg_real frac=Nfaces*EPSILON;
        wpx[fld] = 0.0;
        wpy[fld] = 0.0;
        for(i=0;i<Nfaces;i++){
            gra_x[i] = px[i*Nfield+fld];
            gra_y[i] = py[i*Nfield+fld];
            gra_det[i] = gra_x[i]*gra_x[i] + gra_y[i]*gra_y[i];
        }
        for(i=0;i<Nfaces;i++){
            dg_real w = 1.0;
            for(j=0;j<Nfaces;j++){
                if(i==j){ continue; }
                w = w*gra_det[j];
            }
            frac += w;
            w += EPSILON;
            wpx[fld] += w*gra_x[i];
            wpy[fld] += w*gra_y[i];
#if DEBUG
            printf("fld=%d, f=%d, w=%f\n", fld, i, w);
#endif
        }
        wpx[fld] /= frac;
        wpy[fld] /= frac;
    }
}

void pf_vert_limit(physField *phys){
    const int K = phys->grid->K;
    const int Nfield = phys->Nfield;
    const int nprocs = phys->mesh->nprocs;
    const int Nfaces = phys->cell->Nfaces;
    const int Np = phys->cell->Np;
    const int Nvert = phys->grid->Nv;
    const int Nfp = phys->cell->Nfp;

    parallMesh *mesh = phys->mesh;
    multiReg *region = phys->region;
    dg_real *f_Q = phys->f_Q;

    register int k,n,f,fld;
    /** 1. calculate the cell average value */
    pf_cellMean(phys);
    /* 2. fetch cell info with other processes */
    MPI_Request mpi_out_requests[nprocs];
    MPI_Request mpi_in_requests[nprocs];
    int Nmess;
    /* do sends and recv */
    pf_fetchCellBuffer(phys, mpi_out_requests, mpi_in_requests, &Nmess);

    /* wait for 2. sends and recv buffers */
    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_in_requests, instatus);
    MPI_Waitall(Nmess, mpi_out_requests, instatus);

    /** 2. limited vertex value  */
    int **EToV = phys->grid->EToV;
    /* maximum and minimum of each vertex */
    dg_real *vmax = vector_real_create(Nvert * Nfield);
    dg_real *vmin = vector_real_create(Nvert * Nfield);
    for(n=0;n<Nvert*Nfield;n++){
        vmax[n] = -INFTY;
        vmin[n] =  INFTY;
    }
    /* local vertex */
    for(k=0;k<K;k++){
        for(f=0;f<Nfaces;f++){
            int v = EToV[k][f];
            for(fld=0;fld<Nfield;fld++){
                int sk = v*Nfield+fld;
                vmax[sk] = max(phys->c_Q[k*Nfield+fld], vmax[sk]);
                vmin[sk] = min(phys->c_Q[k*Nfield+fld], vmin[sk]);
            }
        }
    }
    /* adjacent parallel vertex */
    for(n=0;n<mesh->parallCellNum;n++){
        k = mesh->cellIndexIn[n];
        f = mesh->faceIndexIn[n];
        int v1 = EToV[k][f];
        int v2 = EToV[k][(f+1)/Nfaces];
        for(fld=0;fld<Nfield;fld++){
            int sk = v1*Nfield+fld;
            vmax[sk] = max(phys->c_inQ[n*Nfield+fld], vmax[sk]);
            vmin[sk] = min(phys->c_inQ[n*Nfield+fld], vmin[sk]);
            sk = v2*Nfield+fld;
            vmax[sk] = max(phys->c_inQ[n*Nfield+fld], vmax[sk]);
            vmin[sk] = min(phys->c_inQ[n*Nfield+fld], vmin[sk]);
        }
    }
#if DEBUG
    FILE *fp = CreateLog("pf_vert_limit", phys->mesh->procid, phys->mesh->nprocs);
    PrintVector2File(fp, "vmax", vmax, Nvert*Nfield);
    PrintVector2File(fp, "vmin", vmin, Nvert*Nfield);
#endif
    /* limit vertex value of each element */
    for(k=0;k<K;k++){
        for(f=0;f<Nfaces;f++){
            int *fmask = phys->cell->Fmask[f];
            int sk = k*Np*Nfield+ fmask[0]; // node index
            int v = EToV[k][f]; // vertex index
            for(fld=0;fld<Nfield;fld++){
                if(f_Q[sk*Nfield+fld] > vmax[v*Nfield+fld]){
                    f_Q[sk*Nfield+fld] = vmax[v*Nfield+fld];
                }else if (f_Q[sk*Nfield+fld] < vmin[v*Nfield+fld]){
                    f_Q[sk*Nfield+fld] = vmin[v*Nfield+fld];
                }
#if DEBUG
                fprintf(fp, "k=%d, f=%d, fld=%d, f_Q=%f\n", k,f,fld,f_Q[sk*Nfield+fld]);
#endif
            }
        }
    }
    vector_real_free(vmax);
    vector_real_free(vmin);

    dg_real xt[3], yt[3], ft[3];
    dg_real *px = vector_real_create(Nfaces*Nfield);
    dg_real *py = vector_real_create(Nfaces*Nfield);
    for(k=0;k<K;k++){
        double A = 1.0/phys->region->size[k];
        double xc = A * mr_reg_integral(region, k, region->x[k]);
        double yc = A * mr_reg_integral(region, k, region->y[k]);
        xt[0] = xc;
        yt[0] = yc;
        for(f=0;f<Nfaces;f++){
            int *fmask = phys->cell->Fmask[f];
            int v1 = fmask[0]; // local node index
            int v2 = fmask[Nfp-1]; // local node index
            xt[1] = region->x[k][v1];
            yt[1] = region->y[k][v1];
            xt[2] = region->x[k][v2];
            yt[2] = region->y[k][v2];
            v1 = k*Np*Nfield+ fmask[0]; // global node index
            v2 = k*Np*Nfield+ fmask[Nfp-1]; // global node index
            for(fld=0;fld<Nfield;fld++){
                /* 3. gradient of each face */
                ft[0] = phys->c_Q[k*Nfield+fld];
                ft[1] = f_Q[v1*Nfield + fld];
                ft[2] = f_Q[v2*Nfield + fld];
                triGradient(xt, yt, ft, px+f*Nfield+fld, py+f*Nfield+fld);
#if DEBUG
                fprintf(fp, "k=%d, f=%d, fld=%d, px=%f, py=%f\n", k,f,fld,px[f*Nfield+fld],py[f*Nfield+fld] );
#endif
            }
        }
        /* 4. weighted gradient */
        dg_real wpx[Nfield], wpy[Nfield];
        VA_weigrad(Nfaces, Nfield, px, py, wpx, wpy);
#if DEBUG
        for(fld=0;fld<Nfield;fld++){
            fprintf(fp, "k=%d, fld=%d, wpx=%f, wpy=%f\n", k,fld,wpx[fld],wpy[fld]);
        }
#endif
        /* 5. reconstruction */
        for(n=0;n<Np;n++){
            dg_real dx = (dg_real)(region->x[k][n] - xc);
            dg_real dy = (dg_real)(region->y[k][n] - yc);
            for(fld=0;fld<Nfield;fld++){
                dg_real qmean = phys->c_Q[k*Nfield+fld];
                f_Q[n*Nfield+fld] = qmean + dx*wpx[fld] + dy*wpy[fld];
            }
        }
    }
    vector_real_free(px);
    vector_real_free(py);

#if DEBUG
    fclose(fp);
#endif
}