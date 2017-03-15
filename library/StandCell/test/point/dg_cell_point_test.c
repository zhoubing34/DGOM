//
// Created by li12242 on 17/3/15.
//

#include "dg_cell_point_data3.h"
#include "../dg_cell_test_main.h"
#include "dg_cell_point_test.h"

int dg_point_info_test(dg_cell *cell, int verbose){
    int fail = 0;
    if(verbose){
        FILE *fp = fopen("dg_point_info_test.txt", "ws");
        fprintf(fp, "N = %d\n", dg_cell_N(cell));
        fprintf(fp, "Nv = %d\n", dg_cell_Nv(cell));
        fprintf(fp, "Nfaces = %d\n", dg_cell_Nfaces(cell));
        fprintf(fp, "type = %d\n", dg_cell_celltype(cell));
        //print_int_vector2file(fp, "face_type", dg_cell_facetype(cell), dg_cell_Nfaces(cell));
        print_double_vector2file(fp, "vr", dg_cell_vr(cell), dg_cell_Nv(cell));
        print_int_matrix2file(fp, "FToV", dg_cell_FToV(cell), dg_cell_Nfaces(cell), 2);
        fclose(fp);
    }
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_point_node_test(dg_cell *cell, int verbose){
    int fail = 0;
    extern double point_r[NP];
    fail= vector_double_test("dg_point_nood_r_test", dg_cell_r(cell), point_r, dg_cell_Np(cell));
    if(verbose){
        FILE *fp = fopen("dg_point_node_test.txt","ws");
        print_double_vector2file(fp, "r", dg_cell_r(cell), dg_cell_Np(cell));
        fclose(fp);
    }
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_point_vand_matrix_test(dg_cell *cell, int verbose){
    int fail = 0;
    extern double point_V[NP][NP];
    double **V_ext = matrix_double_create(NP, NP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            V_ext[i][j] = point_V[i][j];
        }
    }
    fail = matrix_double_test(__FUNCTION__, dg_cell_V(cell), V_ext, NP, NP);

    if(verbose){
        FILE *fp = fopen("dg_point_vand_matrix_test.txt", "ws");
        print_double_matrix2file(fp, "V", dg_cell_V(cell), NP, NP);
        fclose(fp);
    }
    matrix_double_free(V_ext);
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_point_mass_matrix_test(dg_cell *cell, int verbose){
    int fail = 0;
    extern double point_M[NP][NP];
    double **M_ext = matrix_double_create(NP, NP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            M_ext[i][j] = point_M[i][j];
        }
    }
    fail = matrix_double_test(__FUNCTION__, dg_cell_M(cell), M_ext, NP, NP);

    if(verbose){
        FILE *fp = fopen("dg_point_mass_matrix_test.txt", "ws");
        print_double_matrix2file(fp, "M", dg_cell_M(cell), NP, NP);
        fclose(fp);
    }
    matrix_double_free(M_ext);
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_point_deri_matrix_test(dg_cell *cell, int verbose){
    int fail =0;
    extern double point_Dr[NP][NP];
    double **Dr_ext = matrix_double_create(NP, NP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            Dr_ext[i][j] = point_Dr[i][j];
        }
    }
    fail = matrix_double_test(__FUNCTION__, dg_cell_Dr(cell), Dr_ext, NP, NP);
    if(verbose){
        FILE *fp = fopen("dg_point_deri_matrix_test.txt", "ws");
        print_double_matrix2file(fp, "Dr", dg_cell_Dr(cell), NP, NP);
        fclose(fp);
    }
    matrix_double_free(Dr_ext);
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}


int dg_point_Fmask_test(dg_cell *cell, int verbose){
    int fail = 0;
    const int Nfaces = dg_cell_Nfaces(cell);
    if(verbose){
        FILE *fp = fopen("dg_point_Fmask_test.txt","ws");
        fprintf(fp, "Fmsk = \n");
        int f,n;
        for(f=0;f<Nfaces;f++){
            const int Nfp = dg_cell_Nfp(cell)[f];
            for(n=0;n<Nfp;n++){
                fprintf(fp, "%d ", dg_cell_Fmask(cell)[f][n]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}


int dg_point_LIFT_test(dg_cell *quad, int verbose){
    int fail = 0;
    if(verbose){
        FILE *fp = fopen("dg_point_LIFT_test.txt", "ws");
        print_double_matrix2file(fp, "LIFT", dg_cell_LIFT(quad), NP, NFP);
        fclose(fp);
    }
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}