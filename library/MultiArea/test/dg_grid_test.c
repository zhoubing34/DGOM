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
    int procid,nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    if(verbose){
        /* gen log filename */
        FILE *fp = create_log(__FUNCTION__, procid, nprocs);
        print_int_vector2file(fp, "K", &(grid->K), 1);
        fclose(fp);
    }
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
    int procid,nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    if(verbose){
        /* gen log filename */
        FILE *fp = create_log(__FUNCTION__, procid, nprocs);
        print_int_matrix2file(fp, "EToV", grid->EToV, grid->K, dg_cell_Nv(grid->cell));
        print_int_vector2file(fp, "EToR", grid->EToR, grid->K);
        fclose(fp);
    }
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
    int procid,nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    if(verbose){
        /* gen log filename */
        FILE *fp = create_log(__FUNCTION__, procid, nprocs);
        print_double_vector2file(fp, "vx", grid->vx, grid->Nv);
        print_double_vector2file(fp, "vy", grid->vy, grid->Nv);
        fclose(fp);
    }
    if(!procid) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_grid_connect_test(dg_grid *grid, int verbose){
    int fail = 0;
    int procid,nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    const int Nv = dg_cell_Nv(grid->cell);
    const int K = dg_grid_K(grid);
    if(verbose){
        /* gen log filename */
        FILE *fp = create_log(__FUNCTION__, procid, nprocs);
        print_int_matrix2file(fp, "EToE", grid->EToE, K, Nv);
        print_int_matrix2file(fp, "EToF", grid->EToF, K, Nv);
        print_int_matrix2file(fp, "EToP", grid->EToP, K, Nv);
        fclose(fp);
    }
    if(!procid) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_grid_EToBS_test(dg_grid *grid, int verbose){
    int fail = 0;
    int procid,nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    const int Nv = dg_cell_Nv(grid->cell);
    const int K = dg_grid_K(grid);
    if(verbose){
        /* gen log filename */
        FILE *fp = create_log(__FUNCTION__, procid, nprocs);
        print_int_matrix2file(fp, "EToBS", grid->EToBS, K, Nv);
        fclose(fp);
    }
    if(!procid) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}
