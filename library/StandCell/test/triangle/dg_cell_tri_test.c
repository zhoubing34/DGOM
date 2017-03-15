//
// Created by li12242 on 16/12/16.
//

#include "dg_cell_tri_test.h"
#include "sc_tri_data3.h"

int dg_tri_info_test(dg_cell *cell, int verbose){
    int fail = 0;
    if(verbose){
        FILE *fp = fopen("dg_tri_info_test.txt", "ws");
        fprintf(fp, "N = %d\n", dg_cell_N(cell));
        fprintf(fp, "Nv = %d\n", dg_cell_Nv(cell));
        fprintf(fp, "Nfaces = %d\n", dg_cell_Nfaces(cell));
        fprintf(fp, "type = %d\n", dg_cell_celltype(cell));
        //print_int_vector2file(fp, "face_type", dg_cell_facetype(cell), dg_cell_Nfaces(cell));
        print_double_vector2file(fp, "vr", dg_cell_vr(cell), dg_cell_Nv(cell));
        print_double_vector2file(fp, "vs", dg_cell_vs(cell), dg_cell_Nv(cell));
        print_int_matrix2file(fp, "FToV", dg_cell_FToV(cell), dg_cell_Nfaces(cell), 2);
        fclose(fp);
    }
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_tri_nood_test(dg_cell *tri, int verbose){

    int fail=0;
    extern double tri_r[NP];
    extern double tri_s[NP];

    fail= vector_double_test("dg_tri_nood_r_test", dg_cell_r(tri), tri_r, dg_cell_Np(tri));
    fail= vector_double_test("dg_tri_nood_s_test", dg_cell_s(tri), tri_s, dg_cell_Np(tri));

    if(verbose){
        FILE *fp = fopen("dg_tri_nood_test.txt","ws");
        print_double_vector2file(fp, "r", dg_cell_r(tri), dg_cell_Np(tri));
        print_double_vector2file(fp, "s", dg_cell_s(tri), dg_cell_Np(tri));
        fclose(fp);
    }
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_tri_vand_matrix_test(dg_cell *tri, int verbose){
    int fail = 0;
    extern double tri_V[NP][NP];
    double **V_ext = matrix_double_create(NP, NP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            V_ext[i][j] = tri_V[i][j];
        }
    }
    fail = matrix_double_test(__FUNCTION__, dg_cell_V(tri), V_ext, NP, NP);

    if(verbose){
        FILE *fp = fopen("dg_tri_vand_matrix_test.txt", "ws");
        print_double_matrix2file(fp, "V", dg_cell_V(tri), NP, NP);
        fclose(fp);
    }
    matrix_double_free(V_ext);
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_tri_deri_matrix_test(dg_cell *tri, int verbose){
    int fail =0;
    extern double tri_Dr[NP][NP];
    extern double tri_Ds[NP][NP];
    double **Dr_ext = matrix_double_create(NP, NP);
    double **Ds_ext = matrix_double_create(NP, NP);

    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            Dr_ext[i][j] = tri_Dr[i][j];
            Ds_ext[i][j] = tri_Ds[i][j];
        }
    }

    fail = matrix_double_test(__FUNCTION__, dg_cell_Dr(tri), Dr_ext, NP, NP);
    fail = matrix_double_test(__FUNCTION__, dg_cell_Ds(tri), Ds_ext, NP, NP);

    if(verbose){
        FILE *fp = fopen("dg_tri_deri_matrix_test.txt", "ws");
        print_double_matrix2file(fp, "Dr", dg_cell_Dr(tri), NP, NP);
        print_double_matrix2file(fp, "Ds", dg_cell_Ds(tri), NP, NP);
        fclose(fp);
    }

    matrix_double_free(Dr_ext);
    matrix_double_free(Ds_ext);
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_tri_mass_matrix_test(dg_cell *tri, int verbose){
    int fail = 0;
    extern double tri_M[NP][NP];
    double **M_ext = matrix_double_create(NP, NP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            M_ext[i][j] = tri_M[i][j];
        }
    }
    fail = matrix_double_test(__FUNCTION__, dg_cell_M(tri), M_ext, NP, NP);

    if(verbose){
        FILE *fp = fopen("dg_tri_mass_matrix_test.txt", "ws");
        print_double_matrix2file(fp, "M", dg_cell_M(tri), NP, NP);
        fclose(fp);
    }
    matrix_double_free(M_ext);
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_tri_Fmask_test(dg_cell *cell, int verbose){
    int fail = 0;
    const int Nfaces = dg_cell_Nfaces(cell);
    if(verbose){
        FILE *fp = fopen("dg_tri_Fmask_test.txt","ws");
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
    return fail;
}

int dg_tri_LIFT_test(dg_cell *tri, int verbose){
    int fail = 0;
    extern double tri_LIFT[NP][NFP];
    double **LIFT_ext = matrix_double_create(NP, NFP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NFP;j++){
            LIFT_ext[i][j] = tri_LIFT[i][j];
        }
    }
    fail = matrix_double_test(__FUNCTION__, dg_cell_LIFT(tri), LIFT_ext, NP, NFP);

    if(verbose){
        FILE *fp = fopen("dg_tri_LIFT_test.txt", "ws");
        print_double_matrix2file(fp, "LIFT", dg_cell_LIFT(tri), NP, NFP);
        fclose(fp);
    }

    matrix_double_free(LIFT_ext);
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

//int sc_triVertProj_test(dg_cell *tri, int verbose){
//    int fail = 0;
//    const int Np = dg_cell_Np(tri);
//    extern double tri_VX[NV];
//    extern double tri_VY[NV];
//
//    double x[Np], y[Np];
//
//    dg_cell_proj_vert2node(tri, tri_VX, x);
//    fail = vector_double_test("sc_triVertProj_x_test", x, dg_cell_r(tri), Np);
//
//    dg_cell_proj_vert2node(tri, tri_VY, y);
//    fail = vector_double_test("sc_triVertProj_y_test", y, dg_cell_s(tri), Np);
//    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
//    return fail;
//}