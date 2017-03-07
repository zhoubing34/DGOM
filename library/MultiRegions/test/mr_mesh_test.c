//
// Created by li12242 on 12/21/16.
//
#include "mr_mesh_test.h"

int mr_mesh_connet_test(dg_mesh *mesh, int verbose){
    int fail = 0;
    const int K = mesh->grid->K;
    const int Nfaces = mesh->cell->Nfaces;
    const int Nfp = mesh->cell->Nfp;
    if(verbose){
        FILE *fp = create_log(__FUNCTION__, mesh->procid, mesh->nprocs);
        print_int_matrix2file(fp, "mesh->EToE", mesh->EToE, K, Nfaces);
        print_int_matrix2file(fp, "mesh->EToF", mesh->EToF, K, Nfaces);
        print_int_vector2file(fp, "mesh->vmapM", mesh->vmapM, K * Nfaces * Nfp);
        print_int_vector2file(fp, "mesh->vmapP", mesh->vmapP, K * Nfaces * Nfp);
        fclose(fp);
    }
    const int procid = mesh->procid;
    if(!procid) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int mr_mesh_parallel_test(dg_mesh *mesh, int verbose){
    int fail = 0;
    const int K = mesh->grid->K;
    const int Nfaces = mesh->cell->Nfaces;
    if(verbose){
        FILE *fp = create_log(__FUNCTION__, mesh->procid, mesh->nprocs);
        print_int_matrix2file(fp, "mesh->EToP", mesh->EToP, K, Nfaces);
        print_int_vector2file(fp, "mesh->Npar", mesh->Npar, mesh->nprocs);
        fprintf(fp, "mesh->parallCellNum: %d\n", mesh->parallCellNum);
        print_int_vector2file(fp, "mesh->cellIndexIn", mesh->cellIndexIn, mesh->parallCellNum);
        print_int_vector2file(fp, "mesh->faceIndexIn", mesh->faceIndexIn, mesh->parallCellNum);
        print_int_vector2file(fp, "mesh->cellIndexOut", mesh->cellIndexOut, mesh->parallCellNum);
        fprintf(fp, "mesh->parallNodeNum: %d\n", mesh->parallNodeNum);
        print_int_vector2file(fp, "mesh->nodeIndexOut", mesh->nodeIndexOut, mesh->parallNodeNum);

        fclose(fp);
    }
    const int procid = mesh->procid;
    if(!procid) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int mr_mesh_bc_test(dg_mesh *mesh, int verbose){
    int fail = 0;
    const int K = mesh->grid->K;
    const int Nfaces = mesh->cell->Nfaces;
    if(verbose){
        FILE *fp = create_log(__FUNCTION__, mesh->procid, mesh->nprocs);
        print_int_matrix2file(fp, "mesh->EToBS", mesh->EToBS, K, Nfaces);
        fprintf(fp, "mesh->Nbc: %d\n", mesh->Nbc);
        print_int_vector2file(fp, "mesh->bcind", mesh->bcind, mesh->Nbc);
        fprintf(fp, "mesh->Nobc: %d\n", mesh->Nobc);
        print_int_vector2file(fp, "mesh->obcind", mesh->obcind, mesh->Nobc);
        fclose(fp);
    }
    const int procid = mesh->procid;
    if(!procid) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}
