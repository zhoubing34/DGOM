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
int mr_reg_node_test(multiReg *reg, int verbose){
    int fail = 0;
    const int K = reg->grid->K;
    const int Np = reg->grid->cell->Np;
    if(verbose){
        FILE *fp = create_log(__FUNCTION__, reg->procid, reg->nprocs);
        PrintMatrix2File(fp, "region->x", reg->x, K, Np);
        PrintMatrix2File(fp, "region->y", reg->y, K, Np);
        fclose(fp);
    }
    return fail;
}
/**
 * @brief print the volume factors
 * @param reg
 * @param verbose
 * @return
 */
int mr_reg_volume_factor_test(multiReg *reg, int verbose){
    int fail = 0;
    const int K = reg->grid->K;
    const int Np = reg->grid->cell->Np;
    if(verbose){
        FILE *fp = create_log(__FUNCTION__, reg->procid, reg->nprocs);
        PrintMatrix2File(fp, "region->drdx", reg->drdx, K, Np);
        PrintMatrix2File(fp, "region->drdy", reg->drdy, K, Np);
        PrintMatrix2File(fp, "region->dsdx", reg->dsdx, K, Np);
        PrintMatrix2File(fp, "region->dsdy", reg->dsdy, K, Np);
        PrintMatrix2File(fp, "region->J", reg->J, K, Np);
        fclose(fp);
    }
    return fail;
}
/**
 * @brief print the face factors
 * @param reg
 * @param verbose
 * @return
 */
int mr_reg_face_factor_test(multiReg *reg, int verbose){
    int fail = 0;
    const int K = reg->grid->K;
    const int Nfaces = reg->grid->cell->Nfaces;
    if(verbose){
        FILE *fp = create_log(__FUNCTION__, reg->procid, reg->nprocs);
        PrintMatrix2File(fp, "region->nx", reg->nx, K, Nfaces);
        PrintMatrix2File(fp, "region->ny", reg->ny, K, Nfaces);
        PrintMatrix2File(fp, "region->sJ", reg->sJ, K, Nfaces);
        fclose(fp);
    }
    return fail;
}
/**
 * @brief print the length scale of regions
 * @param reg
 * @param verbose
 * @return
 */
int mr_reg_scale_test(multiReg *reg, int verbose){
    int fail = 0;
    const int K = reg->grid->K;
    if(verbose){
        FILE *fp = create_log(__FUNCTION__, reg->procid, reg->nprocs);
        PrintVector2File(fp, "region->len", reg->len, K);
        PrintVector2File(fp, "region->size", reg->size, K);
        fclose(fp);
    }
    return fail;
}