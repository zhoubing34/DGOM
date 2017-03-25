//
// Created by li12242 on 12/19/16.
//

#include "dg_reg_test.h"

/**
 * @brief print the node coordinate
 * @param reg
 * @param verbose
 * @return
 */
int dg_region_node_test(dg_region *reg, int verbose){
    int fail = 0;
    const int K = dg_grid_K(reg->grid);
    const int Np = dg_cell_Np(reg->cell);
    if(verbose){
        FILE *fp = create_log(__FUNCTION__, reg->procid, reg->nprocs);
        print_double_matrix2file(fp, "region->x", reg->x, K, Np);
        print_double_matrix2file(fp, "region->y", reg->y, K, Np);
        fclose(fp);
    }
    const int procid = reg->procid;
    if(!procid) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}
/**
 * @brief print the volume factors
 * @param reg
 * @param verbose
 * @return
 */
int dg_region_volume_factor_test(dg_region *reg, int verbose){
    int fail = 0;
    const int K = dg_grid_K(reg->grid);
    const int Np = dg_cell_Np(reg->cell);
    if(verbose){
        FILE *fp = create_log(__FUNCTION__, reg->procid, reg->nprocs);
        print_double_matrix2file(fp, "region->drdx", reg->drdx, K, Np);
        print_double_matrix2file(fp, "region->drdy", reg->drdy, K, Np);
        print_double_matrix2file(fp, "region->dsdx", reg->dsdx, K, Np);
        print_double_matrix2file(fp, "region->dsdy", reg->dsdy, K, Np);
        print_double_matrix2file(fp, "region->J", reg->J, K, Np);
        fclose(fp);
    }
    const int procid = reg->procid;
    if(!procid) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}
/**
 * @brief print the face factors
 * @param reg
 * @param verbose
 * @return
 */
int dg_region_face_factor_test(dg_region *reg, int verbose){
    int fail = 0;
    const int K = dg_grid_K(reg->grid);
    const int Nfaces = dg_cell_Nfaces(reg->cell);
    if(verbose){
        FILE *fp = create_log(__FUNCTION__, reg->procid, reg->nprocs);
        print_double_matrix2file(fp, "region->nx", reg->nx, K, Nfaces);
        print_double_matrix2file(fp, "region->ny", reg->ny, K, Nfaces);
        print_double_matrix2file(fp, "region->sJ", reg->sJ, K, Nfaces);
        fclose(fp);
    }
    const int procid = reg->procid;
    if(!procid) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}
/**
 * @brief print the length scale of regions
 * @param reg
 * @param verbose
 * @return
 */
int dg_region_scale_test(dg_region *reg, int verbose){
    int fail = 0;
    const int K = dg_grid_K(reg->grid);
    const int Nfaces = dg_cell_Nfaces(reg->cell);
    if(verbose){
        FILE *fp = create_log(__FUNCTION__, reg->procid, reg->nprocs);
        print_double_vector2file(fp, "region->len", reg->len, K);
        print_double_vector2file(fp, "region->size", reg->size, K);
        print_double_matrix2file(fp, "region->face_size", reg->face_size, K, Nfaces);
        fclose(fp);
    }
    const int procid = reg->procid;
    if(!procid) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_region_face_integral_test(dg_region *region, int verbose){
    int fail = 0;
    const int K = dg_grid_K(region->grid);
    const int Nfaces = dg_cell_Nfaces(region->cell);
    const int Nfield = 2;
    const int Np = dg_cell_Np(region->cell);
    int k,f,n,fld,sk=0;
    dg_real *f_Q = vector_real_create(Nfield*K*Np);
    dg_real *face_Q = vector_real_create(Nfield*K*Nfaces);
    dg_real *face_extQ = vector_real_create(Nfield*K*Nfaces);
    // assignment node value
    for(k=0;k<K;k++){
        for(n=0;n<Np;n++){
            f_Q[sk++] = dg_region_x(region)[k][n];
            f_Q[sk++] = dg_region_y(region)[k][n];
        }
    }
    // exact value
    sk = 0;
    for(k=0;k<K;k++){
        for(f=0;f<Nfaces;f++){
            const int Nfp = dg_cell_Nfp(region->cell)[f];
            int *fmask = dg_cell_Fmask(region->cell)[f];
            double x1 = dg_region_x(region)[k][fmask[0]];
            double x2 = dg_region_x(region)[k][fmask[Nfp-1]];
            face_extQ[sk++] = (x1+x2)*0.5;
            double y1 = dg_region_y(region)[k][fmask[0]];
            double y2 = dg_region_y(region)[k][fmask[Nfp-1]];
            face_extQ[sk++] = (y1+y2)*0.5;
        }
    }
    dg_real ones[Np], face_len[Nfaces];
    for(n=0;n<Np;n++){
        ones[n] = 1.0;
    }

    for(k=0;k<K;k++){
        region->face_integral(region, Nfield, k, f_Q+k*Np*Nfield, face_Q+k*Nfaces*Nfield);
        region->face_integral(region, 1, k, ones, face_len);
        for(f=0;f<Nfaces;f++){
            for(fld=0;fld<Nfield;fld++){
                face_Q[k*Nfaces*Nfield + f*Nfield + fld] /= face_len[f];
            }
        }
    }

    fail = vector_double_test(__FUNCTION__, face_Q, face_extQ, Nfield*K*Nfaces);

    if(verbose){
        FILE *fp = create_log(__FUNCTION__, region->procid, region->nprocs);
        print_double_vector2file(fp, "face_Q", face_Q, Nfield*K*Nfaces);
        print_double_vector2file(fp, "face_extQ", face_extQ, Nfield*K*Nfaces);
        fclose(fp);
    }
    vector_real_free(f_Q);
    vector_real_free(face_Q);
    vector_real_free(face_extQ);
    return fail;
}