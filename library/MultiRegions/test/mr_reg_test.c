//
// Created by li12242 on 12/19/16.
//

#include "mr_reg_test.h"

/**
 * @brief print the node coordinate
 * @param reg
 * @param verbose
 * @return
 */
int mr_reg_node_test(dg_region *reg, int verbose){
    int fail = 0;
    const int K = reg->grid->K;
    const int Np = reg->grid->cell->Np;
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
int mr_reg_volume_factor_test(dg_region *reg, int verbose){
    int fail = 0;
    const int K = reg->grid->K;
    const int Np = reg->grid->cell->Np;
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
int mr_reg_face_factor_test(dg_region *reg, int verbose){
    int fail = 0;
    const int K = reg->grid->K;
    const int Nfaces = reg->grid->cell->Nfaces;
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
int mr_reg_scale_test(dg_region *reg, int verbose){
    int fail = 0;
    const int K = reg->grid->K;
    if(verbose){
        FILE *fp = create_log(__FUNCTION__, reg->procid, reg->nprocs);
        print_double_vector2file(fp, "region->len", reg->len, K);
        print_double_vector2file(fp, "region->size", reg->size, K);
        fclose(fp);
    }
    const int procid = reg->procid;
    if(!procid) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}