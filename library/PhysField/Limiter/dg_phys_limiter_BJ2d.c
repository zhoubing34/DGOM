//
// Created by li12242 on 17/2/28.
//

#include "dg_phys_limiter_BJ2d.h"
#define DEBUG 0
#if DEBUG
#include "Utility/unit_test.h"
#endif

/**
 * @brief integral averaged values of each faces
 * @param[in] phys_info physical field structure.
 * @param[in,out] f_mean integral averaged values of each faces
 */
static void dg_phys_local_face_mean(dg_phys_info *phys_info, dg_real *f_mean){

    dg_cell *cell = phys_info->cell;
    dg_region *region = phys_info->region;
    const int Nfield = phys_info->Nfield;
    const int Nfaces = dg_cell_Nfaces(cell);
    const int Np = dg_cell_Np(cell);
    const int K = dg_grid_K(phys_info->grid);

    dg_real *f_Q = phys_info->f_Q;

    register int k,f,fld,sk=0;
    for(k=0;k<K;k++){
        region->face_integral(region, Nfield, k, f_Q+k*Np*Nfield, f_mean+k*Nfaces*Nfield);
        for(f=0;f<Nfaces;f++){
            double ds = region->face_size[k][f];
            for(fld=0;fld<Nfield;fld++)
                f_mean[sk++] /= ds;
        }
    }
    return;
}

static void dg_phys_weight_face_mean(dg_phys_info *phys_info, dg_real *f_mean){

    dg_region *region = phys_info->region;
    dg_mesh *mesh = phys_info->mesh;
    dg_grid *grid = phys_info->grid;

    const int Nfield = phys_info->Nfield;
    const int Nfaces = dg_cell_Nfaces(phys_info->cell);
    const int K = dg_grid_K(grid);
    const int procid = dg_grid_procid(grid);
    const int nprocs = dg_grid_nprocs(grid);

    register int k,f,n,fld;
    /* calculate the center coordinate */
    dg_real *xc = vector_real_create(K);
    dg_real *yc = vector_real_create(K);
    for(k=0;k<K;k++){
        const double Area = 1.0/region->size[k];
        region->vol_integral(region, 1, k, region->x[k], xc+k);
        region->vol_integral(region, 1, k, region->y[k], yc+k);
        xc[k] *= Area;
        yc[k] *= Area;
    }
    const int Nfetchfaces = dg_mesh_NfetchFace(mesh);
    dg_real *xc_in = vector_real_create(Nfetchfaces);
    dg_real *yc_in = vector_real_create(Nfetchfaces);

    MPI_Request xc_out_requests[nprocs];
    MPI_Request xc_in_requests[nprocs];
    MPI_Request yc_out_requests[nprocs];
    MPI_Request yc_in_requests[nprocs];

    int Nmess;
    Nmess = mesh->fetch_cell_buffer(mesh, 1, xc, xc_in, xc_out_requests, xc_in_requests);
    Nmess = mesh->fetch_cell_buffer(mesh, 1, yc, yc_in, yc_out_requests, yc_in_requests);

    dg_real *c_Q = phys_info->c_Q; ///< cell mean value;
    int **EToE = dg_grid_EToE(grid);
    int **EToP = dg_grid_EToP(grid);
    double **x = region->x;
    double **y = region->y;
    for(k=0;k<K;k++){
        double xf[Nfaces], yf[Nfaces];
        region->face_integral(region, 1, k, x[k], xf);
        region->face_integral(region, 1, k, y[k], yf);
        for(f=0;f<Nfaces;f++){

            int e = EToE[k][f];
            int p = EToP[k][f];
            if(p!=procid) { continue; }
            double ds = region->face_size[k][f];
            xf[f] /= ds;
            yf[f] /= ds;

            double d1 = (xf[f] - xc[k])*(xf[f] - xc[k]) + (yf[f] - yc[k])*(yf[f] - yc[k]);
            double d2 = (xf[f] - xc[e])*(xf[f] - xc[e]) + (yf[f] - yc[e])*(yf[f] - yc[e]);
            dg_real w1 = d2/(d1 + d2);
            dg_real w2 = d1/(d1 + d2);

            int sk = k*Nfaces*Nfield + f*Nfield;
            for(fld=0;fld<Nfield;fld++)
                f_mean[sk+fld] = (c_Q[k*Nfield+fld]*w1 + c_Q[e*Nfield+fld]*w2);
        }
    }
    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, xc_in_requests, instatus);
    MPI_Waitall(Nmess, xc_out_requests, instatus);
    MPI_Waitall(Nmess, yc_in_requests, instatus);
    MPI_Waitall(Nmess, yc_out_requests, instatus);

    for(n=0;n<Nfetchfaces;n++){
        k = mesh->CBFToK[n];
        f = mesh->CBFToF[n];

        double ds = region->face_size[k][f];
        dg_real xc_next = xc_in[n];
        dg_real yc_next = yc_in[n];

        double xf[Nfaces], yf[Nfaces], zf[Nfaces];
        region->face_integral(region, 1, k, x[k], xf);
        region->face_integral(region, 1, k, y[k], yf);
        xf[f] /= ds;
        yf[f] /= ds;

        double d1 = (xf[f] - xc[k])*(xf[f] - xc[k]) + (yf[f] - yc[k])*(yf[f] - yc[k]);
        double d2 = (xf[f] - xc_next)*(xf[f] - xc_next) + (yf[f] - yc_next)*(yf[f] - yc_next);
        dg_real w1 = d2/(d1 + d2);
        dg_real w2 = d1/(d1 + d2);

        int sk = k*Nfaces*Nfield + f*Nfield;
        for(fld=0;fld<Nfield;fld++){
            f_mean[sk+fld] = (c_Q[k*Nfield+fld]*w1 + phys_info->c_recvQ[n*Nfield+fld]*w2);
        }
    }

    vector_real_free(xc_in);
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
static void dg_phys_gradient(dg_phys_info *phys, dg_real *px, dg_real *py){
    dg_region *region = phys->region;
    const int K = dg_grid_K(phys->grid);
    const int Nfield = phys->Nfield;
    const int Nfaces = dg_cell_Nfaces(phys->cell);

    register int k,f,fld,sk;
    dg_real *f_mean = vector_real_create(K*Nfaces*Nfield);
    dg_phys_local_face_mean(phys, f_mean);

    for(k=0;k<K;k++){
        double A = 1.0/dg_region_size(region)[k];

        for(fld=0;fld<Nfield;fld++){
            sk = k*Nfield+fld;
            px[sk] = 0;
            py[sk] = 0;
        }

        for(f=0;f<Nfaces;f++){
            const double ds = dg_region_face_size(region)[k][f];
            const double nx = dg_region_nx(region)[k][f];
            const double ny = dg_region_ny(region)[k][f];

            for(fld=0;fld<Nfield;fld++){
                sk = k*Nfield+fld;
                int sf = k*Nfaces*Nfield + f*Nfield;
                px[sk] += f_mean[sf+fld]*nx*ds*A;
                py[sk] += f_mean[sf+fld]*ny*ds*A;
            }
        }
    }
    vector_real_free(f_mean);
    return;
}

/**
 *
 * @param phys
 * @param cell_max
 * @param cell_min
 */
static void dg_phys_adjacent_cellinfo(dg_phys_info *phys, dg_real *cell_max, dg_real *cell_min){
    dg_grid *grid = phys->grid;
    const int K = dg_grid_K(grid);
    const int Nfield = phys->Nfield;
    const int procid = dg_grid_procid(grid);
    const int nprocs = dg_grid_nprocs(grid);
    const int Nfaces = dg_cell_Nfaces(phys->cell);

    register int k,n,f,fld;

    MPI_Request mpi_out_requests[nprocs];
    MPI_Request mpi_in_requests[nprocs];
    int Nmess;

    /* do sends and recv */
    Nmess = phys->fetch_cell_buffer(phys, mpi_out_requests, mpi_in_requests);

    /* local cell loop */
    int **EToE = dg_grid_EToE(grid);
    int **EToP = dg_grid_EToP(grid);
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
    dg_mesh *mesh = phys->mesh;
    const int NfetchFace = dg_mesh_NfetchFace(mesh);
    for(n=0;n<NfetchFace;n++){
        k = mesh->CBFToK[n];
        for(fld=0;fld<Nfield;fld++){
            dg_real c_next = phys->c_recvQ[n*Nfield+fld];
            int sk = k*Nfield+fld;
            cell_max[sk] = max(cell_max[sk], c_next);
            cell_min[sk] = min(cell_min[sk], c_next);
        }
    }
}

/**
 * @brief obtain Barth-Jeperson slope limiter.
 * @param phys_info
 * @param Cmax
 * @param Cmin
 * @param beta
 * @param psi
 */
static void dg_phys_BJ_limiter(dg_phys_info *phys_info,
                               dg_real *Cmax, dg_real *Cmin,
                               double beta, dg_real *psi){

    const int K = phys_info->grid->K;
    const int Nfield = phys_info->Nfield;
    const int Np = dg_cell_Np(phys_info->cell);

    register int k,fld,n;
    for(k=0;k<K;k++){
        dg_real *f_Q = phys_info->f_Q + k*Np*Nfield; // variable of k-th cell
        for(fld=0;fld<Nfield;fld++){
            int sk = k*Nfield+fld;
            dg_real qmax = Cmax[sk];
            dg_real qmin = Cmin[sk];
            dg_real qmean = phys_info->c_Q[sk];
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
 * The slope limiter will act on each physical variables and reconstruct the scalar distribution.
 * @param[in,out] phys pointer to dg_phys structure;
 * @param[in] tind trouble cell indicator (0 for smooth; 1 for trouble cell)
 * @param[in] beta
 */
void dg_phys_limiter_BJ2d(dg_phys_info *phys, int *tind, double beta){

    const int K = dg_grid_K(phys->grid);
    const int Nfield = phys->Nfield;
    const int Np = dg_cell_Np(phys->cell);
    register int k,n,fld;

#if DEBUG
    FILE *fp = create_log("dg_phys_limiter_BJ2d", phys->mesh->procid, phys->mesh->nprocs);
#endif

    /* 1. calculate the cell average value */
    phys->cell_mean(phys);

    /* 2. fetch cell info with other processes */
    dg_real cell_max[K*Nfield], cell_min[K*Nfield];
    dg_phys_adjacent_cellinfo(phys, cell_max, cell_min);
#if DEBUG
    print_double_vector2file(fp, "cell_max", cell_max, K*Nfield);
    print_double_vector2file(fp, "cell_min", cell_min, K*Nfield);
#endif
    /* 3. calculate the unlimited gradient */
    dg_real px[K*Nfield], py[K*Nfield];
    dg_phys_gradient(phys, px, py);

#if DEBUG
    print_double_vector2file(fp, "px", px, K*Nfield);
    print_double_vector2file(fp, "py", py, K*Nfield);
#endif

    /* 4. calculate the limited results */
    dg_real psi[K*Nfield];
    dg_phys_BJ_limiter(phys, cell_max, cell_min, beta, psi);

#if DEBUG
    print_int_vector2file(fp, "tind", tind, K*Nfield);
#endif

    /* 6. reconstruction */
    dg_region *region = phys->region;
    for(k=0;k<K;k++){
        double xc, yc;
        const double Area = 1.0/region->size[k];
        region->vol_integral(region, 1, k, region->x[k], &xc);
        region->vol_integral(region, 1, k, region->y[k], &yc);
        xc *= Area;
        yc *= Area;

        dg_real *f_Q = phys->f_Q + k*Np*Nfield; // variable of k-th cell
        for(fld=0;fld<Nfield;fld++){
            if( tind[k*Nfield+fld] == 0 ) continue;
            int sk = k*Nfield+fld;
            /* compute the limited gradient */
            dg_real qx = px[sk] * psi[sk];
            dg_real qy = py[sk] * psi[sk];

            /* reconstruction of each element */
            dg_real qmean = phys->c_Q[sk];
            for(n=0;n<Np;n++){
                dg_real dx = (dg_real)(dg_region_x(region)[k][n] - xc);
                dg_real dy = (dg_real)(dg_region_y(region)[k][n] - yc);

                f_Q[n*Nfield+fld] = qmean + dx*qx + dy*qy;
            }
        }
    }

#if DEBUG
    fclose(fp);
#endif
    return;
}