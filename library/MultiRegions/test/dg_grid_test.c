//
// Created by li12242 on 12/18/16.
//

#include "dg_grid_test.h"

/**
 * @brief print the number of element in each process
 * @param grid
 * @param verbose
 * @return
 */
int dg_grid_Ktol_test(dg_grid *grid, int verbose){
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
int dg_grid_EToV_test(dg_grid *grid, int verbose){
    int fail = 0;
    if(verbose){
        /* gen log filename */
        FILE *fp = create_log(__FUNCTION__, grid->procid, grid->nprocs);
        print_int_matrix2file(fp, "EToV", grid->EToV, grid->K, grid->cell->Nv);
        print_int_vector2file(fp, "EToR", grid->EToR, grid->K);
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
int dg_grid_vertex_test(dg_grid *grid, int verbose){
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

int dg_grid_connect_test(dg_grid *grid, int verbose){
    int fail = 0;
    const int Nv = dg_cell_Nv(grid->cell);
    const int K = dg_grid_K(grid);
    if(verbose){
        /* gen log filename */
        FILE *fp = create_log(__FUNCTION__, grid->procid, grid->nprocs);
        print_int_matrix2file(fp, "EToE", grid->EToE, K, Nv);
        print_int_matrix2file(fp, "EToF", grid->EToF, K, Nv);
        print_int_matrix2file(fp, "EToP", grid->EToP, K, Nv);
        fclose(fp);
    }
    const int procid = grid->procid;
    if(!procid) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_grid_EToBS_test(dg_grid *grid, int verbose){
    int fail = 0;
    const int Nv = dg_cell_Nv(grid->cell);
    const int K = dg_grid_K(grid);
    if(verbose){
        /* gen log filename */
        FILE *fp = create_log(__FUNCTION__, grid->procid, grid->nprocs);
        print_int_matrix2file(fp, "EToBS", grid->EToBS, K, Nv);
        fclose(fp);
    }
    return fail;
}
