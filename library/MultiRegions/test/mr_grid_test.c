//
// Created by li12242 on 12/18/16.
//

#include "mr_grid_test.h"

/**
 * @brief print the number of element in each process
 * @param grid
 * @param verbose
 * @return
 */
int mr_grid_Ktol_test(geoGrid *grid, int verbose){
    int fail = 0;
    if(verbose){
        /* gen log filename */
        FILE *fp = create_log(__FUNCTION__, grid->procid, grid->nprocs);
        print_int_vector2file(fp, "K", &(grid->K), 1);
        fclose(fp);
    }
    const int procid = grid->procid;
    if(!procid) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

/**
 * @brief print the EToV of each process
 * @param grid
 * @param verbose
 * @return
 */
int mr_grid_EToV_test(geoGrid *grid, int verbose){
    int fail = 0;
    if(verbose){
        /* gen log filename */
        FILE *fp = create_log(__FUNCTION__, grid->procid, grid->nprocs);
        print_int_matrix2file(fp, "EToV", grid->EToV, grid->K, grid->cell->Nv);
        fclose(fp);
    }
    const int procid = grid->procid;
    if(!procid) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}
/**
 * @brief print the vertex coordinate
 * @param grid
 * @param verbose
 * @return
 */
int mr_grid_vertex_test(geoGrid *grid, int verbose){
    int fail = 0;
    if(verbose){
        /* gen log filename */
        FILE *fp = create_log(__FUNCTION__, grid->procid, grid->nprocs);
        print_double_vector2file(fp, "vx", grid->vx, grid->Nv);
        print_double_vector2file(fp, "vy", grid->vy, grid->Nv);
        fclose(fp);
    }
    const int procid = grid->procid;
    if(!procid) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}
