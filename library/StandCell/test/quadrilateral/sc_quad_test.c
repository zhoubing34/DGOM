//
// Created by li12242 on 16/12/17.
//

#include "sc_quad_test.h"
#include "sc_quad_data3.h"

int sc_quadCoor_test(stdCell *quad, int verbose){
    int fail = 0;
    const int Np = quad->Np;
    extern double quad_r[NP];
    extern double quad_s[NP];

    fail=Vector_test("sc_quadCoor_r_test", quad->r, quad_r, Np, 0);
    fail=Vector_test("sc_quadCoor_s_test", quad->s, quad_s, Np, 0);

    if(verbose){
        FILE *fp = fopen("sc_quadCoor_test.txt","w");
        PrintVector2File(fp, "r", quad->r, Np);
        PrintVector2File(fp, "s", quad->s, Np);
        fclose(fp);
    }
    return fail;
}

int sc_quadVandMatrix_test(stdCell *quad, int verbose){
    int fail = 0;
    extern double quad_V[NP][NP];
    double **V_ext = Matrix_create(NP, NP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            V_ext[i][j] = quad_V[i][j];
        }
    }
    fail = Matrix_test("sc_quadVandMatrix_test", quad->V, V_ext, NP, NP, 0);

    if(verbose){
        FILE *fp = fopen("sc_quadVandMatrix_test.txt", "w");
        PrintMatrix2File(fp, "V", quad->V, NP, NP);
        fclose(fp);
    }

    Matrix_free(V_ext);
    return fail;
}

int sc_quadMassMatrix_test(stdCell *quad, int verbose){
    int fail = 0;
    extern double quad_M[NP][NP];
    double **M_ext = Matrix_create(NP, NP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            M_ext[i][j] = quad_M[i][j];
        }
    }
    fail = Matrix_test("sc_quadMassMatrix_test", quad->M, M_ext, NP, NP, 0);

    if(verbose){
        FILE *fp = fopen("sc_quadMassMatrix_test.txt", "w");
        PrintMatrix2File(fp, "M", quad->M, NP, NP);
        fclose(fp);
    }
    Matrix_free(M_ext);
    return fail;
}

int sc_quadDeriMatrix_test(stdCell *quad, int verbose){
    int fail =0;
    extern double quad_Dr[NP][NP];
    extern double quad_Ds[NP][NP];
    double **Dr_ext = Matrix_create(NP, NP);
    double **Ds_ext = Matrix_create(NP, NP);

    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            Dr_ext[i][j] = quad_Dr[i][j];
            Ds_ext[i][j] = quad_Ds[i][j];
        }
    }

    fail = Matrix_test("sc_quadDr_test", quad->Dr, Dr_ext, NP, NP, 0);
    fail = Matrix_test("sc_quadDs_test", quad->Dr, Dr_ext, NP, NP, 0);

    if(verbose){
        FILE *fp = fopen("sc_quadDeriMatrix_test.txt", "w");
        PrintMatrix2File(fp, "Dr", quad->Dr, NP, NP);
        PrintMatrix2File(fp, "Ds", quad->Ds, NP, NP);
        fclose(fp);
    }

    Matrix_free(Dr_ext);
    Matrix_free(Ds_ext);
    return fail;
}


int sc_quadLIFT_test(stdCell *quad, int verbose){
    int fail = 0;
    extern double quad_LIFT[NP][NFP];
    double **LIFT_ext = Matrix_create(NP, NFP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NFP;j++){
            LIFT_ext[i][j] = quad_LIFT[i][j];
        }
    }

    fail = Matrix_test("sc_quadLIFT_test", quad->LIFT, LIFT_ext, NP, NFP, 0);

    if(verbose){
        FILE *fp = fopen("sc_quadLIFT_test.txt", "w");
        PrintMatrix2File(fp, "LIFT", quad->LIFT, NP, NFP);
        fclose(fp);
    }

    Matrix_free(LIFT_ext);
    return fail;
}