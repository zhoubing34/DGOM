//
// Created by li12242 on 16/12/16.
//

#include "sc_tri_test.h"
#include "sc_tri_data3.h"
#include "StandCell/test/sc_test.h"

int sc_triCoor_test(stdCell *tri, int verbose){

    int fail=0;
    extern double tri_r[NP];
    extern double tri_s[NP];

    fail=Vector_test("sc_triCoor_r_test", tri->r, tri_r, tri->Np, 0);
    fail=Vector_test("sc_triCoor_s_test", tri->s, tri_s, tri->Np, 0);

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
    double **V_ext = Matrix_create(NP, NP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            V_ext[i][j] = tri_V[i][j];
        }
    }
    fail = Matrix_test("sc_triVandMatrix_test", tri->V, V_ext, NP, NP, 0);

    if(verbose){
        FILE *fp = fopen("sc_triVandMatrix_test.txt", "w");
        PrintMatrix2File(fp, "V", tri->V, NP, NP);
        fclose(fp);
    }

    Matrix_free(V_ext);
    return fail;
}

int sc_triMassMatrix_test(stdCell *tri, int verbose){
    int fail = 0;
    extern double tri_M[NP][NP];
    double **M_ext = Matrix_create(NP, NP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            M_ext[i][j] = tri_M[i][j];
        }
    }
    fail = Matrix_test("sc_triMassMatrix_test", tri->M, M_ext, NP, NP, 0);

    if(verbose){
        FILE *fp = fopen("sc_triMassMatrix_test.txt", "w");
        PrintMatrix2File(fp, "M", tri->M, NP, NP);
        fclose(fp);
    }
    Matrix_free(M_ext);
    return fail;
}

int sc_triDeriMatrix_test(stdCell *tri, int verbose){
    int fail =0;
    extern double tri_Dr[NP][NP];
    extern double tri_Ds[NP][NP];
    double **Dr_ext = Matrix_create(NP, NP);
    double **Ds_ext = Matrix_create(NP, NP);

    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            Dr_ext[i][j] = tri_Dr[i][j];
            Ds_ext[i][j] = tri_Ds[i][j];
        }
    }

    fail = Matrix_test("sc_triDr_test", tri->Dr, Dr_ext, NP, NP, 0);
    fail = Matrix_test("sc_triDs_test", tri->Dr, Dr_ext, NP, NP, 0);

    if(verbose){
        FILE *fp = fopen("sc_triDeriMatrix_test.txt", "w");
        PrintMatrix2File(fp, "Dr", tri->Dr, NP, NP);
        PrintMatrix2File(fp, "Ds", tri->Ds, NP, NP);
        fclose(fp);
    }

    Matrix_free(Dr_ext);
    Matrix_free(Ds_ext);
    return fail;
}

int sc_triLIFT_test(stdCell *tri, int verbose){
    int fail = 0;
    extern double tri_LIFT[NP][NFP];
    double **LIFT_ext = Matrix_create(NP, NFP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NFP;j++){
            LIFT_ext[i][j] = tri_LIFT[i][j];
        }
    }

    fail = Matrix_test("sc_triLIFT_test", tri->LIFT, LIFT_ext, NP, NFP, 0);

    if(verbose){
        FILE *fp = fopen("sc_triLIFT_test.txt", "w");
        PrintMatrix2File(fp, "LIFT", tri->LIFT, NP, NFP);
        fclose(fp);
    }

    Matrix_free(LIFT_ext);
    return fail;
}