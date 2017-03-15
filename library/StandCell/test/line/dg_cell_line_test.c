//
// Created by li12242 on 17/3/15.
//

#include "dg_cell_line_data3.h"
#include "../dg_cell_test_main.h"
#include "dg_cell_line_test.h"

int dg_line_info_test(dg_cell *cell, int verbose){
    int fail = 0;
    if(verbose){
        FILE *fp = fopen("dg_line_info_test.txt", "ws");
        fprintf(fp, "N = %d\n", dg_cell_N(cell));
        fprintf(fp, "Nv = %d\n", dg_cell_Nv(cell));
        fprintf(fp, "Nfaces = %d\n", dg_cell_Nfaces(cell));
        fprintf(fp, "type = %d\n", dg_cell_celltype(cell));
        //print_int_vector2file(fp, "face_type", dg_cell_facetype(cell), dg_cell_Nfaces(cell));
        print_double_vector2file(fp, "vr", dg_cell_vr(cell), dg_cell_Nv(cell));
        print_int_matrix2file(fp, "FToV", dg_cell_FToV(cell), dg_cell_Nfaces(cell), 1);
        fclose(fp);
    }
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_line_node_test(dg_cell *cell, const int verbose){
    int fail = 0;
    extern double line_r[NP];
    fail= vector_double_test("dg_tri_nood_r_test", dg_cell_r(cell), line_r, dg_cell_Np(cell));
    if(verbose){
        FILE *fp = fopen("dg_line_nood_test.txt","ws");
        print_double_vector2file(fp, "r", dg_cell_r(cell), dg_cell_Np(cell));
        fclose(fp);
    }
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_line_vand_matrix_test(dg_cell *cell, int verbose){
    int fail = 0;
    extern double line_V[NP][NP];
    double **V_ext = matrix_double_create(NP, NP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            V_ext[i][j] = line_V[i][j];
        }
    }
    fail = matrix_double_test(__FUNCTION__, dg_cell_V(cell), V_ext, NP, NP);

    if(verbose){
        FILE *fp = fopen("dg_line_vand_matrix_test.txt", "ws");
        print_double_matrix2file(fp, "V", dg_cell_V(cell), NP, NP);
        fclose(fp);
    }
    matrix_double_free(V_ext);
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_line_mass_matrix_test(dg_cell *cell, int verbose){
    int fail = 0;
    extern double line_M[NP][NP];
    double **M_ext = matrix_double_create(NP, NP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            M_ext[i][j] = line_M[i][j];
        }
    }
    fail = matrix_double_test(__FUNCTION__, dg_cell_M(cell), M_ext, NP, NP);

    if(verbose){
        FILE *fp = fopen("dg_line_vand_matrix_test.txt", "ws");
        print_double_matrix2file(fp, "V", dg_cell_M(cell), NP, NP);
        fclose(fp);
    }
    matrix_double_free(M_ext);
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_line_deri_matrix_test(dg_cell *cell, int verbose){
    int fail =0;
    extern double line_Dr[NP][NP];
    double **Dr_ext = matrix_double_create(NP, NP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            Dr_ext[i][j] = line_Dr[i][j];
        }
    }
    fail = matrix_double_test(__FUNCTION__, dg_cell_Dr(cell), Dr_ext, NP, NP);
    if(verbose){
        FILE *fp = fopen("dg_line_deri_matrix_test.txt", "ws");
        print_double_matrix2file(fp, "Dr", dg_cell_Dr(cell), NP, NP);
        fclose(fp);
    }
    matrix_double_free(Dr_ext);
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_line_Fmask_test(dg_cell *cell, int verbose){
    int fail = 0;
    const int Nfaces = dg_cell_Nfaces(cell);
    if(verbose){
        FILE *fp = fopen("dg_line_Fmask_test.txt","ws");
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


int dg_line_LIFT_test(dg_cell *quad, int verbose){
    int fail = 0;
    extern double line_LIFT[NP][NFP];
    double **LIFT_ext = matrix_double_create(NP, NFP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NFP;j++){
            LIFT_ext[i][j] = line_LIFT[i][j];
        }
    }

    fail = matrix_double_test(__FUNCTION__, dg_cell_LIFT(quad), LIFT_ext, NP, NFP);

    if(verbose){
        FILE *fp = fopen("dg_line_LIFT_test.txt", "ws");
        print_double_matrix2file(fp, "LIFT", dg_cell_LIFT(quad), NP, NFP);
        fclose(fp);
    }

    matrix_double_free(LIFT_ext);
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}