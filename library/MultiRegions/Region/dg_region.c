#include "dg_region.h"
#include "dg_region_volumInfo.h"
#include "dg_region_surfInfo.h"

#define DEBUG 0

static void dg_region_node(dg_region *region);
static void dg_region_scale2d(dg_region *region);
static void dg_region_scale3d(dg_region *region);
static void dg_region_vol_integral(dg_region *region, int Nfield, int k, dg_real *f_Q,
                                   dg_real *c_Q);
static void dg_region_face_integral(dg_region *region, int Nfield, int k, dg_real *f_Q,
                                    dg_real *face_Q);
/**
 * @brief functions of creating structure dg_region
 */
typedef struct dg_region_creator{
    void (*set_nood)(dg_region *reg);
    void (*set_volumInfo)(dg_region *reg);
    void (*set_surfInfo)(dg_region *reg);
    void (*set_scale)(dg_region *reg);
}dg_region_creator;

static const dg_region_creator region2d_creator={
        dg_region_node,
        dg_reg_volumInfo2d,
        dg_reg_surfInfo2d,
        dg_region_scale2d,
};

static const dg_region_creator region3d_creator={
        dg_region_node,
        dg_reg_volumInfo3d,
        dg_reg_surfInfo3d,
        dg_region_scale3d,
};

dg_region* dg_region_create(dg_grid *grid){

    dg_region *region = (dg_region *)calloc(1, sizeof(dg_region));
    /* basic infomation */
    region->procid = grid->procid;
    region->nprocs = grid->nprocs;
    region->cell = grid->cell;
    region->grid = grid;

    const dg_region_creator *creator;
    switch (dg_cell_celltype(grid->cell)){
        case TRIANGLE:
            creator = &region2d_creator; break;
        case QUADRIL:
            creator = &region2d_creator; break;
        default:
            fprintf(stderr, "%s (%d): Unknown cell type %d\n",
                    __FUNCTION__, __LINE__, dg_cell_celltype(grid->cell));
            exit(-1);
    }
    region->vol_integral = dg_region_vol_integral;
    region->face_integral = dg_region_face_integral;
    creator->set_nood(region);
    creator->set_volumInfo(region);
    creator->set_surfInfo(region);
    creator->set_scale(region);

    return region;
}

void dg_region_free(dg_region *region){
    matrix_double_free(region->x);
    matrix_double_free(region->y);
    matrix_double_free(region->z);
    matrix_double_free(region->J);
    /* volume geometry */
    matrix_double_free(region->drdx);
    matrix_double_free(region->drdy);
    matrix_double_free(region->dsdx);
    matrix_double_free(region->dsdy);

    matrix_double_free(region->nx);
    matrix_double_free(region->ny);
    matrix_double_free(region->sJ);

    /* volume/area size and length len */
    vector_double_free(region->size);
    vector_double_free(region->len);

    free(region);
    return;
}

/**
 * @brief create node coordinate for 2d (triangle and quadrilateral)
 * @param[in,out] region multi-regions object
 */
static void dg_region_node(dg_region *region){
    dg_cell *cell = region->cell;
    dg_grid *grid = region->grid;
    const int Np = dg_cell_Np(cell);
    const int Nv = dg_cell_Nv(cell);
    const int K = dg_grid_K(grid);
    int **EToV = dg_grid_EToV(grid);
    double *vx = dg_grid_vx(grid);
    double *vy = dg_grid_vy(grid);
    double *vz = dg_grid_vz(grid);

    region->x = matrix_double_create(K, Np);
    region->y = matrix_double_create(K, Np);
    region->z = matrix_double_create(K, Np);

    int k,i;
    double gx[Nv], gy[Nv], gz[Nv];
    for(k=0;k<K;k++){
        for(i=0;i<Nv;i++){
            gx[i] = vx[EToV[k][i]];
            gy[i] = vy[EToV[k][i]];
            gz[i] = vz[EToV[k][i]];
        }
        cell->proj_vert2node(cell, 1, gx, region->x[k]);
        cell->proj_vert2node(cell, 1, gy, region->y[k]);
        cell->proj_vert2node(cell, 1, gz, region->z[k]);
    }
    return;
}

static void dg_region_scale2d(dg_region *region){
    const int Np = dg_cell_Np(region->cell);
    const int K = dg_grid_K(region->grid);
    const int Nfaces = dg_cell_Nfaces(region->cell);

    region->size = vector_double_create(K); // volume or area
    region->len = vector_double_create(K); // length of element
    region->face_size = matrix_double_create(K, Nfaces);
    int k,i;
    /* initialize ones */
    dg_real ones[Np];
    for(i=0;i<Np;i++) {ones[i] = 1.0;}

    // elemental size
    for(k=0;k<K;k++){
        dg_real area;
        region->vol_integral(region, 1, k, ones, &area);
        region->size[k] = area;
        region->len[k] = sqrt(area/M_PI);
        region->face_integral(region, 1, k, ones, region->face_size[k]);
    }
    return;
}

static void dg_region_scale3d(dg_region *region){
    const int Np = dg_cell_Np(region->cell);
    const int K = dg_grid_K(region->grid);

    region->size = vector_double_create(K); // volume or area
    region->len = vector_double_create(K); // length of element
    int k,i;
    /* initialize ones */
    double ones[Np];
    for(i=0;i<Np;i++) {ones[i] = 1.0;}

    // elemental size
    for(k=0;k<K;k++){
        double area;
        region->vol_integral(region, 1, k, ones, &area);
        region->size[k] = area;
        region->len[k] = pow(area*3.0/4./M_PI, 1.0/3);
    }
    return;
}

/**
 * @brief calculate the integral value in kth cell.
 *
 * @param[in] region pointer to a dg_region structure;
 * @param[in] Nfield number of physical field;
 * @param[in] k index of integral element;
 * @param[in] f_Q value on interpolation points;
 * @param[out] c_Q integral value of each physical field;
 */
static void dg_region_vol_integral(dg_region *region, int Nfield, int k,
                                   dg_real *f_Q, dg_real *c_Q){

    const double *J = region->J[k];
    const double *w = dg_cell_w(region->cell);
    const int Np = dg_cell_Np(region->cell);

    register int i,fld;
    // initialize
    for(fld=0;fld<Nfield;fld++){
        c_Q[fld] = 0;
    }
    // vol_integral at all nodes
    for(i=0;i<Np;i++){
        for(fld=0;fld<Nfield;fld++){
            c_Q[fld] += J[i]*w[i]*f_Q[i*Nfield+fld];
        }
    }
    return;
}
/**
 * @brief calculate the integral value at each face in kth cell.
 * @param[in] region pointer to a dg_region structure;
 * @param[in] Nfield number of physical field;
 * @param[in] k index of integral element;
 * @param[in] f_Q value on interpolation points;
 * @param[out] face_Q integral value on each face;
 */
static void dg_region_face_integral(dg_region *region, int Nfield, int k,
                                    dg_real *f_Q, dg_real *face_Q){
    dg_cell *cell = region->cell;
    const int Nfaces = dg_cell_Nfaces(cell);

    register int f,n,fld;
    // initialize
    for(fld=0;fld<Nfield*Nfaces;fld++){
        face_Q[fld] = 0;
    }
    // vol_integral at all faces
    for(f=0;f<Nfaces;f++){
        const int Nfp = dg_cell_Nfp(cell)[f];
        int *fmask = dg_cell_Fmask(cell)[f];
        const double sJ = region->sJ[k][f];
        for(n=0;n<Nfp;n++){
            const double ws = dg_cell_ws(cell)[f][n];
            for(fld=0;fld<Nfield;fld++){
                face_Q[f*Nfield + fld] += sJ*ws*f_Q[fmask[n]*Nfield+fld];
            }
        }
    }
}