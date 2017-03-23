//
// Created by li12242 on 16/12/17.
//

#include "dg_cell_quad_test.h"
#include "sc_quad_data3.h"

int dg_quad_info_test(dg_cell *cell, int verbose){
    int fail = 0;
    if(verbose){
        FILE *fp = fopen("dg_quad_info_test.txt", "ws");
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

int dg_quad_nood_test(dg_cell *quad, int verbose){
    int fail = 0;
    const int Np = dg_cell_Np(quad);
    extern double quad_r[NP];
    extern double quad_s[NP];

    fail= vector_double_test(__FUNCTION__, dg_cell_r(quad), quad_r, Np);
    fail= vector_double_test(__FUNCTION__, dg_cell_s(quad), quad_s, Np);

    if(verbose){
        FILE *fp = fopen("dg_quad_nood_test.txt","ws");
        print_double_vector2file(fp, "r", dg_cell_r(quad), Np);
        print_double_vector2file(fp, "s", dg_cell_s(quad), Np);
        fclose(fp);
    }
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_quad_vand_matrix_test(dg_cell *quad, int verbose){
    int fail = 0;
    extern double quad_V[NP][NP];
    double **V_ext = matrix_double_create(NP, NP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            V_ext[i][j] = quad_V[i][j];
        }
    }
    fail = matrix_double_test(__FUNCTION__, dg_cell_V(quad), V_ext, NP, NP);

    if(verbose){
        FILE *fp = fopen("dg_quad_vand_matrix_test.txt", "ws");
        print_double_matrix2file(fp, "V", dg_cell_V(quad), NP, NP);
        fclose(fp);
    }

    matrix_double_free(V_ext);
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_quad_deri_matrix_test(dg_cell *quad, int verbose){
    int fail =0;
    extern double quad_Dr[NP][NP];
    extern double quad_Ds[NP][NP];
    double **Dr_ext = matrix_double_create(NP, NP);
    double **Ds_ext = matrix_double_create(NP, NP);

    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            Dr_ext[i][j] = quad_Dr[i][j];
            Ds_ext[i][j] = quad_Ds[i][j];
        }
    }

    fail = matrix_double_test(__FUNCTION__, dg_cell_Dr(quad), Dr_ext, NP, NP);
    fail = matrix_double_test(__FUNCTION__, dg_cell_Ds(quad), Ds_ext, NP, NP);

    if(verbose){
        FILE *fp = fopen("dg_quad_deri_matrix_test.txt", "ws");
        print_double_matrix2file(fp, "Dr", dg_cell_Dr(quad), NP, NP);
        print_double_matrix2file(fp, "Ds", dg_cell_Ds(quad), NP, NP);
        fclose(fp);
    }

    matrix_double_free(Dr_ext);
    matrix_double_free(Ds_ext);
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_quad_mass_matrix_test(dg_cell *quad, int verbose){
    int fail = 0;
    extern double quad_M[NP][NP];
    double **M_ext = matrix_double_create(NP, NP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            M_ext[i][j] = quad_M[i][j];
        }
    }
    fail = matrix_double_test(__FUNCTION__, dg_cell_M(quad), M_ext, NP, NP);

    if(verbose){
        FILE *fp = fopen("dg_quad_mass_matrix_test.txt", "ws");
        print_double_matrix2file(fp, "M", dg_cell_M(quad), NP, NP);
        fclose(fp);
    }
    matrix_double_free(M_ext);
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_quad_Fmask_test(dg_cell *cell, int verbose){
    int fail = 0;
    const int Nfaces = dg_cell_Nfaces(cell);
    if(verbose){
        FILE *fp = fopen("dg_quad_Fmask_test.txt","ws");
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


int dg_quad_LIFT_test(dg_cell *quad, int verbose){
    int fail = 0;
    extern double quad_LIFT[NP][NFP];
    double **LIFT_ext = matrix_double_create(NP, NFP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NFP;j++){
            LIFT_ext[i][j] = quad_LIFT[i][j];
        }
    }

    fail = matrix_double_test("dg_quad_LIFT_test", dg_cell_LIFT(quad), LIFT_ext, NP, NFP);

    if(verbose){
        FILE *fp = fopen("dg_quad_LIFT_test.txt", "ws");
        print_double_matrix2file(fp, "LIFT", dg_cell_LIFT(quad), NP, NFP);
        fclose(fp);
    }

    matrix_double_free(LIFT_ext);
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}


int dg_quad_vert_proj_test(dg_cell *quad, int verbose){
    int fail=0;
    extern double quad_VX[NV];
    extern double quad_VY[NV];
    const int Np = dg_cell_Np(quad);
    double x[Np], y[Np];
    dg_cell_proj_vert2node(quad, 1, quad_VX, x);
    fail = vector_double_test("dg_quad_vert_proj_test", x, dg_cell_r(quad), Np);

    dg_cell_proj_vert2node(quad, 1, quad_VY, y);
    fail = vector_double_test("dg_quad_vert_proj_test", y, dg_cell_s(quad), Np);
    if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}