//
// Created by li12242 on 16/12/17.
//

#include "sc_quad_test.h"
#include "sc_quad_data3.h"

int sc_quadCoor_test(dg_cell *quad, int verbose){
    int fail = 0;
    const int Np = quad->Np;
    extern double quad_r[NP];
    extern double quad_s[NP];

    fail= vector_double_test("sc_quadCoor_r_test", quad->r, quad_r, Np);
    fail= vector_double_test("sc_quadCoor_s_test", quad->s, quad_s, Np);

    if(verbose){
        FILE *fp = fopen("sc_quadCoor_test.txt","w");
        print_double_vector2file(fp, "r", quad->r, Np);
        print_double_vector2file(fp, "s", quad->s, Np);
        fclose(fp);
    }
    return fail;
}

int sc_quadFmask_test(dg_cell *cell, int verbose){
    int fail = 0;
    const int Nfaces = dg_cell_Nfaces(cell);
    if(verbose){
        FILE *fp = fopen("sc_quadFmask_test.txt","w");
        fprintf(fp, "Fmsk = \n");
        int f,n;
        for(f=0;f<Nfaces;f++){
            for(n=0;n<dg_cell_Nfp(cell, f);n++){
                fprintf(fp, "%d ", cell->Fmask[f][n]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
    return fail;
}

int sc_quadVandMatrix_test(dg_cell *quad, int verbose){
    int fail = 0;
    extern double quad_V[NP][NP];
    double **V_ext = matrix_double_create(NP, NP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            V_ext[i][j] = quad_V[i][j];
        }
    }
    fail = matrix_double_test("sc_quadVandMatrix_test", quad->V, V_ext, NP, NP);

    if(verbose){
        FILE *fp = fopen("sc_quadVandMatrix_test.txt", "w");
        print_double_matrix2file(fp, "V", quad->V, NP, NP);
        fclose(fp);
    }

    matrix_double_free(V_ext);
    return fail;
}

int sc_quadMassMatrix_test(dg_cell *quad, int verbose){
    int fail = 0;
    extern double quad_M[NP][NP];
    double **M_ext = matrix_double_create(NP, NP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NP;j++){
            M_ext[i][j] = quad_M[i][j];
        }
    }
    fail = matrix_double_test("sc_quadMassMatrix_test", quad->M, M_ext, NP, NP);

    if(verbose){
        FILE *fp = fopen("sc_quadMassMatrix_test.txt", "w");
        print_double_matrix2file(fp, "M", quad->M, NP, NP);
        fclose(fp);
    }
    matrix_double_free(M_ext);
    return fail;
}

int sc_quadDeriMatrix_test(dg_cell *quad, int verbose){
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

    fail = matrix_double_test("sc_quadDr_test", quad->Dr, Dr_ext, NP, NP);
    fail = matrix_double_test("sc_quadDs_test", quad->Ds, Ds_ext, NP, NP);

    if(verbose){
        FILE *fp = fopen("sc_quadDeriMatrix_test.txt", "w");
        print_double_matrix2file(fp, "Dr", quad->Dr, NP, NP);
        print_double_matrix2file(fp, "Ds", quad->Ds, NP, NP);
        fclose(fp);
    }

    matrix_double_free(Dr_ext);
    matrix_double_free(Ds_ext);
    return fail;
}


int sc_quadLIFT_test(dg_cell *quad, int verbose){
    int fail = 0;
    extern double quad_LIFT[NP][NFP];
    double **LIFT_ext = matrix_double_create(NP, NFP);
    int i,j;
    for(i=0;i<NP;i++){
        for(j=0;j<NFP;j++){
            LIFT_ext[i][j] = quad_LIFT[i][j];
        }
    }

    fail = matrix_double_test("sc_quadLIFT_test", quad->LIFT, LIFT_ext, NP, NFP);

    if(verbose){
        FILE *fp = fopen("sc_quadLIFT_test.txt", "w");
        print_double_matrix2file(fp, "LIFT", quad->LIFT, NP, NFP);
        fclose(fp);
    }

    matrix_double_free(LIFT_ext);
    return fail;
}


int sc_quadVertProj_test(dg_cell *quad, int verbose){
    int fail=0;
    extern double quad_VX[NV];
    extern double quad_VY[NV];

    double x[quad->Np], y[quad->Np];
    dg_cell_proj_vert2node(quad, quad_VX, x);
    fail = vector_double_test("sc_quadVertProj_test", x, quad->r, quad->Np);

    dg_cell_proj_vert2node(quad, quad_VY, y);
    fail = vector_double_test("sc_quadVertProj_test", y, quad->s, quad->Np);
    return fail;
}