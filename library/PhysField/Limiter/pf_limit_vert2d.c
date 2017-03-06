//
// Created by li12242 on 17/2/28.
//

#include "pf_limit_vert2d.h"

static void triGradient(dg_real *x, dg_real *y, dg_real *f, dg_real *px, dg_real *py){
    dg_real a[4], b[2];
    a[0] = x[1] - x[0]; a[1] = y[1] - y[0];
    a[2] = x[2] - x[0]; a[3] = y[2] - y[0];
    b[0] = f[1] - f[0]; b[1] = f[2] - f[0];

    dg_real det = a[0]*a[3] - a[1]*a[2];
    *px = ( b[0]*a[3] - b[1]*a[1])/det;
    *py = (-b[0]*a[2] + b[1]*a[0])/det;
    return;
}

static void VA_weigrad(const int Nfaces, const int Nfield,
                       dg_real *px, dg_real *py,
                       dg_real *wpx, dg_real *wpy)
{

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