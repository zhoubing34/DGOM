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
        PrintIntVector2File(fp, "K", &(grid->K), 1);
        fclose(fp);
    }
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
        PrintIntMatrix2File(fp, "EToV", grid->EToV, grid->K, grid->cell->Nv);
        fclose(fp);
    }
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
        PrintVector2File(fp, "vx", grid->vx, grid->cell->Nv);
        PrintVector2File(fp, "vy", grid->vy, grid->cell->Nv);
        fclose(fp);
    }
    return fail;
}
