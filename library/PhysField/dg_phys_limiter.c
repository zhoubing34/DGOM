//
// Created by li12242 on 17/3/23.
//

#include "dg_phys_limiter.h"
#include "PhysField/Limiter/dg_phys_limiter_BJ2d.h"
#include "PhysField/Limiter/dg_phys_indicator_edge.h"

static void dg_phys_set_limiter(dg_phys_limiter *dg_limit, Limiter_Type type);
static void dg_phys_set_indicator(dg_phys_limiter *limiter, Indicator_Type type);
static void dg_phys_limit_trouble_cell(struct dg_phys_limiter *limiter, dg_phys_info *info, double parameter);

dg_phys_limiter* dg_phys_limiter_create(dg_phys_info *info){
    dg_phys_limiter *limiter = (dg_phys_limiter*) calloc(1, sizeof(dg_phys_limiter));

    const int K = dg_grid_K(info->grid);
    const int Nfield = info->Nfield;
    limiter->indicator = vector_int_create(K*Nfield);
    limiter->set_limiter = dg_phys_set_limiter;
    limiter->set_indicator = dg_phys_set_indicator;
    /* default limiter */
    limiter->limiter_func = dg_phys_limiter_BJ2d;
    /* default indicator */
    limiter->indicator_func = dg_phys_indicator_edge;
    /* limiter method */
    limiter->limit_trouble_cell = dg_phys_limit_trouble_cell;
    return limiter;
}

void dg_phys_limiter_free(dg_phys_limiter *limiter){
    vector_int_free(limiter->indicator);
    return;
}

/**
 * @brief limit the result with slope limiter function.
 * @details firstly search for the trouble cell, and store these index into indicator of
 * dg_phys_limiter structure, then limit these cells with specific limiter function.
 *
 * @param limiter pointer to dg_phys_limiter structure;
 * @param info pointer to dg_phys_info structure;
 * @param parameter parameter for the the slope limiter;
 */
static void dg_phys_limit_trouble_cell(struct dg_phys_limiter *limiter, dg_phys_info *info, double parameter){
    limiter->indicator_func(info, limiter->indicator);
    limiter->limiter_func(info, limiter->indicator, parameter);
    return;
}
/**
 * @brief set the indicator function.
 * @param limiter
 * @param type
 */
static void dg_phys_set_indicator(dg_phys_limiter *limiter, Indicator_Type type){
    switch (type){
        case EDGE_INDICATOR:
            limiter->indicator_func = dg_phys_indicator_edge;
            break;
        default:
            fprintf(stderr, "%s (%d): Unknown limiter\n", __FUNCTION__, __LINE__);
            exit(-1);
    }
    return;
}
/**
 * @brief set the slope limiter function.
 * @param limiter
 * @param type
 */
static void dg_phys_set_limiter(dg_phys_limiter *limiter, Limiter_Type type){
    switch (type){
        case BJ_LIMITER:
            limiter->limiter_func = dg_phys_limiter_BJ2d;
            break;
        default:
            fprintf(stderr, "%s (%d): Unknown limiter\n", __FUNCTION__, __LINE__);
            exit(-1);
    }
    return;
}
