//
// Created by li12242 on 16/12/16.
//

#include "sc_tri_test.h"
#include "sc_tri_data3.h"

int sc_triCoor_test(stdCell *tri, int verbose){

    int fail=0;
    extern double tri_r[NP];
    extern double tri_s[NP];

    fail= vector_double_test("sc_triCoor_r_test", tri->r, tri_r, tri->Np);
    fail= vector_double_test("sc_triCoor_s_test", tri->s, tri_s, tri->Np);

    if(verbose){
        FILE *fp = fopen("sc_triCoor_test.txt","w");
        PrintVector2File(fp, "r", tri->r, tri->Np);
        PrintVector2File(fp, "s", tri->s, tri->Np);
        fclose(fp);
    }
    return fail;
}

int sc_Fmask(stdCell *tri, int verbose){
    int fail = 0;
    return fail;
}

int sc_triVandMatrix_test(stdCell *tri, int verbose){
    int fail = 0;
    extern double tri_V[NP][NP];
    double **V_ext = matrix_double_create(NP, NP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            V_ext[i][j] = tri_V[i][j];
        }
    }
    fail = matrix_double_test("sc_triVandMatrix_test", tri->V, V_ext, NP, NP);

    if(verbose){
        FILE *fp = fopen("sc_triVandMatrix_test.txt", "w");
        PrintMatrix2File(fp, "V", tri->V, NP, NP);
        fclose(fp);
    }

    matrix_double_free(V_ext);
    return fail;
}

int sc_triMassMatrix_test(stdCell *tri, int verbose){
    int fail = 0;
    extern double tri_M[NP][NP];
    double **M_ext = matrix_double_create(NP, NP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            M_ext[i][j] = tri_M[i][j];
        }
    }
    fail = matrix_double_test("sc_triMassMatrix_test", tri->M, M_ext, NP, NP);

    if(verbose){
        FILE *fp = fopen("sc_triMassMatrix_test.txt", "w");
        PrintMatrix2File(fp, "M", tri->M, NP, NP);
        fclose(fp);
    }
    matrix_double_free(M_ext);
    return fail;
}

int sc_triDeriMatrix_test(stdCell *tri, int verbose){
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

    fail = matrix_double_test("sc_triDr_test", tri->Dr, Dr_ext, NP, NP);
    fail = matrix_double_test("sc_triDs_test", tri->Ds, Ds_ext, NP, NP);

    if(verbose){
        FILE *fp = fopen("sc_triDeriMatrix_test.txt", "w");
        PrintMatrix2File(fp, "Dr", tri->Dr, NP, NP);
        PrintMatrix2File(fp, "Ds", tri->Ds, NP, NP);
        fclose(fp);
    }

    matrix_double_free(Dr_ext);
    matrix_double_free(Ds_ext);
    return fail;
}

int sc_triLIFT_test(stdCell *tri, int verbose){
    int fail = 0;
    extern double tri_LIFT[NP][NFP];
    double **LIFT_ext = matrix_double_create(NP, NFP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NFP;j++){
            LIFT_ext[i][j] = tri_LIFT[i][j];
        }
    }

    fail = matrix_double_test("sc_triLIFT_test", tri->LIFT, LIFT_ext, NP, NFP);

    if(verbose){
        FILE *fp = fopen("sc_triLIFT_test.txt", "w");
        PrintMatrix2File(fp, "LIFT", tri->LIFT, NP, NFP);
        fclose(fp);
    }

    matrix_double_free(LIFT_ext);
    return fail;
}

int sc_triVertProj_test(stdCell *tri, int verbose){
    int fail = 0;

    extern double tri_VX[NV];
    extern double tri_VY[NV];

    double x[tri->Np], y[tri->Np];

    sc_vertProj(tri, tri_VX, x);
    fail = vector_double_test("sc_triVertProj_x_test", x, tri->r, tri->Np);

    sc_vertProj(tri, tri_VY, y);
    fail = vector_double_test("sc_triVertProj_y_test", y, tri->s, tri->Np);

    return fail;
}