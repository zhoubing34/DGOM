#include <StandCell/sc_stdcell.h>
#include "mr_reg.h"
#include "mr_grid.h"

/* create node coordinate for multiReg object */
static void mr_reg_nodeCoor(multiReg *region);

/* create node coordinate for 2d multiReg object */
static void mr_reg_nodeCoor2d(multiReg *region);

/* create node coordinate for 3d multiReg object */
static void mr_reg_nodeCoor3d(multiReg *region);

/* calculate the volume factors (drdx, dsdx, drdy and dsdy) for regions */
static void mr_reg_volumInfo(multiReg *region);

/* calculate the volume factors for 2d region */
static void mr_reg_volumInfo2d(multiReg *region);

/* calculate the volume factors for 3d region */
static void mr_reg_volumInfo3d(multiReg *region);

/* calculate the volume/area and the length scale of each element */
static void mr_reg_volScale(multiReg *region);

/* calculate the outward normal vector and jacobi coefficient on faces */
static void mr_reg_surfInfo2d(multiReg *region);


/**
 * @brief
 * create of multi-region object.
 *
 * @param [in] grid geometry grid object
 * @return region multi-region object
 *
 * @note
 * postcondition: user should call @ref mr_reg_free to free the memory manually.
 */
multiReg* mr_reg_create(geoGrid *grid){

    /* allocation */
    multiReg *region = (multiReg *)calloc(1, sizeof(multiReg));

    region->dim = grid->dim;
    region->procid = grid->procid;
    region->nprocs = grid->nprocs;

    /* standard element */
    region->cell = grid->cell;
    /* geometry grid */
    region->grid = grid;

    /* node coordinate, x,y and z */
    mr_reg_nodeCoor(region);

    /* volume factor (drdx, dsdx, drdy and dsdy) and J */
    mr_reg_volumInfo(region);

    /* volume/area size and length len */
    mr_reg_volScale(region);

    /* face factors (nx, ny, sJ) */
    mr_reg_surfInfo2d(region);

    return region;
};

/**
 * @brief free the memory of multiReg
 * @param[in] region multi-regions object
 */
void mr_reg_free(multiReg *region){
    const int dim = region->dim;
    /* nodes coordinate */
    switch (dim){
        case 2:
            Matrix_free(region->x);
            Matrix_free(region->y);
            break;
        case 3:
            Matrix_free(region->x);
            Matrix_free(region->y);
            Matrix_free(region->z);
            break;
        default:
            printf("MultiRegions (mr_reg_free): Wrong dimensions %d.\n", dim);
            exit(-1);
    }

    /* Jacobi coefficient */
    Matrix_free(region->J);
    /* volume geometry */
    Matrix_free(region->drdx);
    Matrix_free(region->drdy);
    Matrix_free(region->dsdx);
    Matrix_free(region->dsdy);

    /* surface geometry */
    Matrix_free(region->nx);
    Matrix_free(region->ny);
    Matrix_free(region->sJ);

    /* volume/area size and length len */
    Vector_free(region->size);
    Vector_free(region->len);

    free(region);
}

/**
 * @brief create node coordinate for multiReg object
 * @param[in,out] region multi-regions object
 */
static void mr_reg_nodeCoor(multiReg *region){

    const int dim = region->cell->dim;
    switch (dim){
        case 2:
            mr_reg_nodeCoor2d(region); break;
        case 3:
            mr_reg_nodeCoor3d(region); break;
        default:
            printf("MultiRegions (mr_reg_nodeCoor): Wrong dimensions %d.\n", dim);
            exit(-1);
    }
    return;
}

/**
 * @brief create node coordinate for 2d (triangle and quadrilateral)
 * @param[in,out] region multi-regions object
 */
static void mr_reg_nodeCoor3d(multiReg *region){
    stdCell *cell = region->cell;
    geoGrid *grid = region->grid;
    const int Np = cell->Np;
    const int K = grid->K;
    const int Nv = cell->Nv;
    int **EToV = grid->EToV;
    double *vx = grid->vx;
    double *vy = grid->vy;
    double *vz = grid->vz;

    region->x = Matrix_create(K, Np);
    region->y = Matrix_create(K, Np);
    region->z = Matrix_create(K, Np);

    int k,i;
    double gx[Nv], gy[Nv], gz[Nv];
    for(k=0;k<K;k++){
        for(i=0;i<Nv;i++){
            gx[i] = vx[EToV[k][i]];
            gy[i] = vy[EToV[k][i]];
            gz[i] = vz[EToV[k][i]];
        }
        sc_vertProj(cell, gx, region->x[k]);
        sc_vertProj(cell, gy, region->y[k]);
        sc_vertProj(cell, gz, region->z[k]);
    }
    return;
}

/**
 * @brief create node coordinate for 2d region (triangle and quadrilateral)
 * @param[in,out] region multi-regions object
 */
static void mr_reg_nodeCoor2d(multiReg *region){
    stdCell *cell = region->cell;
    geoGrid *grid = region->grid;
    const int Np = cell->Np;
    const int K = grid->K;
    const int Nv = cell->Nv;
    int **EToV = grid->EToV;
    double *vx = grid->vx;
    double *vy = grid->vy;

    region->x = Matrix_create(K, Np);
    region->y = Matrix_create(K, Np);

    int k,i;
    double gx[Nv], gy[Nv];
    for(k=0;k<K;k++){
        for(i=0;i<Nv;i++){
            gx[i] = vx[EToV[k][i]];
            gy[i] = vy[EToV[k][i]];
        }
        sc_vertProj(cell, gx, region->x[k]);
        sc_vertProj(cell, gy, region->y[k]);
    }
    return;
}


/**
 * @brief  calculate the volume factors (drdx, dsdx, drdy and dsdy) for regions
 * @param[in,out] region multi-region object
 */
static void mr_reg_volumInfo(multiReg *region) {
    const int dim = region->cell->dim;
    switch (dim){
        case 2:
            mr_reg_volumInfo2d(region); break;
        case 3:
            mr_reg_volumInfo3d(region); break;
        default:
            printf("MultiRegions (mr_reg_volumInfo): Wrong dimensions %d.\n", dim);
            exit(-1);
    }
    return;
}

/**
 * @brief calculate the volume factors (drdx, dsdx, drdy and dsdy) for 3d region
 * @param[in,out] region multi-regions object
 * @todo
 */
static void mr_reg_volumInfo3d(multiReg *region){

    return;
}

/**
 * @brief calculate the volume factors (drdx, dsdx, drdy and dsdy) for 2d region (triangle and quadrilateral)
 * @details
 * Get the geometric factor in region, the derived geometric factor on each node is calculated as
 * \f[ \left. \frac{\partial x}{\partial r} \right|_{\mathbf{r}_i} =
 * \left. \frac{\partial \mathbf{\varphi}^T_j}{\partial r} \right|_{\mathbf{r}_i} \mathbf{x}_j \f]
 * where \f$ \mathbf{\varphi} \f$ is the vector of Lagrangian basis functions. The same as
 * \f$ \frac{\partial x}{\partial s}, \frac{\partial y}{\partial r} \frac{\partial y}{\partial s} \f$
 *
 * The Jacobian matrix is used to transfer the differentials of (x,y) to (r,s) with the chain rule:
 * \f[ \begin{bmatrix} dx \cr dy \end{bmatrix} =
 * \begin{bmatrix}
 * \frac{\partial x}{\partial r} & \frac{\partial x}{\partial s} \cr
 * \frac{\partial y}{\partial r} & \frac{\partial y}{\partial s} \cr
 * \end{bmatrix}
 * \begin{bmatrix} dr \cr ds \end{bmatrix} =
 * \mathbf{J}^T \begin{bmatrix} dr \cr ds \end{bmatrix} \f]
 *
 * The corresponding Jacobian and inverse Jacobian is
 * \f[ \mathbf{J} = \frac{\partial (x, y)}{\partial (r, s)} =
 * \begin{bmatrix}
 * \frac{\partial x}{\partial r} & \frac{\partial x}{\partial s} \cr
 * \frac{\partial y}{\partial r} & \frac{\partial y}{\partial s} \cr
 * \end{bmatrix} = \begin{bmatrix} J_{11} & J_{12} \cr
 * J_{21} & J_{22} \end{bmatrix}, \quad
 * \mathbf{J}^{-1} = \frac{\partial (r, s)}{\partial (x,y)} =
 * \begin{bmatrix}
 * \frac{\partial r}{\partial x} & \frac{\partial s}{\partial x} \cr
 * \frac{\partial r}{\partial y} & \frac{\partial s}{\partial y} \cr
 * \end{bmatrix} = \frac{1}{J} \begin{bmatrix} J_{22} & -J_{12} \cr
 * -J_{21} & J_{11} \end{bmatrix}, \quad \f]
 * where \f$ J = det(\mathbf{J}) = J_{11}J_{22} - J_{12}J_{21}\f$
 *
 * @param[in,out] region multi-regions object
 */
static void mr_reg_volumInfo2d(multiReg *region){
    int k,n;
    const int Np = region->cell->Np;
    double **Dr = region->cell->Dr;
    double **Ds = region->cell->Ds;
    const int K = region->grid->K;
    double dxdr[Np], dxds[Np], dydr[Np], dyds[Np];

    // allocation
    region->J = Matrix_create(K, Np);
    region->drdx = Matrix_create(K, Np);
    region->drdy = Matrix_create(K, Np);
    region->dsdx = Matrix_create(K, Np);
    region->dsdy = Matrix_create(K, Np);

    double **drdx = region->drdx;
    double **drdy = region->drdy;
    double **dsdx = region->dsdx;
    double **dsdy = region->dsdy;
//    region->Nvgeo = 4;
//    size_t sz = (size_t)region->Nvgeo*K*Np;
//    region->vgeo = (real*) calloc(sz, sizeof(real));
//    real *vgeo = region->vgeo;

    for(k=0;k<K;k++){

        // calculate the dxdr,dxds,dydr,dyds
        Matrix_multiply(Np, Np, 1, Dr[0], region->x[k], dxdr);
        Matrix_multiply(Np, Np, 1, Ds[0], region->x[k], dxds);
        Matrix_multiply(Np, Np, 1, Dr[0], region->y[k], dydr);
        Matrix_multiply(Np, Np, 1, Ds[0], region->y[k], dyds);

        for(n=0;n<Np;n++){
            double jtemp = -dxds[n]*dydr[n] + dxdr[n]*dyds[n];

            region->J[k][n] = jtemp;
            if(region->J[k][n]<0)
                printf("MultiRegions (mr_reg_volumInfo2d): Warning: J[%d][%d] = %lg\n",
                       k, n, region->J[k][n]);

            /* inverted Jacobian matrix for coordinate mapping */
            drdx[k][n] =  dyds[n]/jtemp; // drdx
            drdy[k][n] = -dxds[n]/jtemp; // drdy
            dsdx[k][n] = -dydr[n]/jtemp; // dsdx
            dsdy[k][n] =  dxdr[n]/jtemp; // dsdy
        }
    }
    return;
}

/**
 * @brief calculate the outward normal vector and jacobi coefficient on faces
 * @param[in] region multi-regions
 */
static void mr_reg_surfInfo2d(multiReg *region){
    int k,f;
    const int K = region->grid->K;
    const int Nfaces = region->cell->Nfaces;
    const int Nfp = region->cell->Nfp;

    double **x = region->x;
    double **y = region->y;
    int **Fmask = region->cell->Fmask;

    double **nx = Matrix_create(K, Nfaces);
    double **ny = Matrix_create(K, Nfaces);
    double **sJ = Matrix_create(K, Nfaces);

    region->nx = nx;
    region->ny = ny;
    region->sJ = sJ;

    for(k=0;k<K;k++){
        for(f=0;f<Nfaces;f++){
            double x1 = x[k][Fmask[f][0]];
            double x2 = x[k][Fmask[f][Nfp-1]];
            double y1 = y[k][Fmask[f][0]];
            double y2 = y[k][Fmask[f][Nfp-1]];

            nx[k][f] =  (y2-y1);
            ny[k][f] = -(x2-x1);

            sJ[k][f] = sqrt(nx[k][f]*nx[k][f]+ny[k][f]*ny[k][f]);
            nx[k][f] /= sJ[k][f];
            ny[k][f] /= sJ[k][f];
            sJ[k][f] /= 2.;
        }
    }

    return;
}

/**
 * @brief calculate the volume/area and the length scale of each element
 * @param[in,out] region multi-region object
 */
static void mr_reg_volScale(multiReg *region){
    stdCell *cell = region->cell;
    geoGrid *grid = region->grid;
    const int Np = cell->Np;
    const int K = grid->K;

    region->size = Vector_create(K); // volume or area
    region->len = Vector_create(K); // length of element
    int k,i;
    /* initialize ones */
    double ones[Np];
    for(i=0;i<Np;i++)
        ones[i] = 1.0;

    // elemental size
    for(k=0;k<K;k++){
        region->size[k] = mr_reg_integral(region, k, ones);
    }

    // elemental scale
    switch (cell->dim){
        case 2: // 2d area
            for(k=0;k<K;k++){
                region->len[k] = sqrt(region->size[k]/M_PI);
            }
            break;
        case 3: // 3d volume
            for(k=0;k<K;k++){
                region->len[k] = pow(region->size[k]*3.0/4./M_PI, 1.0/3);
            }
            break;
        default:
            printf("MultiRegions (mr_reg_volScale): Wrong dimensions %d.\n", cell->dim);
            exit(-1);
    }
}

/**
 * @brief integral in the specific element
 *
 * @param[in] region multi-region object
 * @param[in] ind index of integral element
 * @param[in] nodalVal value on interpolation points
 * @return integral integral value
 */
double mr_reg_integral(multiReg *region, int ind, double *nodalVal){
    double integral = 0;

    stdCell *cell = region->cell;
    double *J = region->J[ind];
    const int Np = cell->Np;

    int i;
    for(i=0;i<Np;i++){
        integral += J[i]*cell->wv[i]*nodalVal[i];
    }
    return integral;
}