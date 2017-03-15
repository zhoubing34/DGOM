//
// Created by li12242 on 17/3/14.
//

#include "dg_cell_face.h"

void dg_cell_face_free(dg_cell_face *face){
    vector_int_free(face->Nfp);
    array_double_free(face->ws);
    array_int_free(face->Fmask);
    matrix_double_free(face->LIFT);
    free(face);
    return;
}

/**
 * @brief
 * set the properties of Nfptotal and Nfp in dg_cell_face structure.
 * @param cell
 * @param face
 */
static void dg_cell_face_info(dg_cell *cell, dg_cell_face *face){
    const int N = dg_cell_N(cell);
    const int Nfaces = dg_cell_Nfaces(cell);

    dg_cell_type *face_type = dg_cell_facetype(cell);
    int *Nfp = vector_int_create(Nfaces);
    int f,Nfptotal=0;
    for(f=0;f<Nfaces;f++){
        dg_cell *face_cell = dg_cell_creat(N, face_type[f]);
        Nfp[f] = dg_cell_Np(face_cell); // node number on fth face
        Nfptotal += Nfp[f];
        dg_cell_free(face_cell);
    }
    face->Nfptotal = Nfptotal;
    face->Nfp = Nfp;
    return;
}

static double **dg_cell_face_LIFT(dg_cell *cell, int Nfptotal, int **Fmask){
    const int N = dg_cell_N(cell);
    const int Np = dg_cell_Np(cell);
    const int Nfaces = dg_cell_Nfaces(cell);

    double **LIFT = matrix_double_create(Np, Nfptotal);
    dg_cell_type *face_type = dg_cell_facetype(cell);
    double *Mes = vector_double_create(Np*Nfptotal);
    int f,n,m,sk=0;
    for(f=0;f<Nfaces;f++){
        dg_cell *face_cell = dg_cell_creat(N, face_type[f]);
        int Nface_node = dg_cell_Np(face_cell);
        /* assignment of Mes */
        for(n=0;n<Nface_node;n++){     /* column index of M */
            for(m=0;m<Nface_node;m++){ /* row index of M */
                int row = Fmask[f][m];   /* row index of Mes */
                Mes[row*Nfptotal + sk] = dg_cell_M(face_cell)[m][n];
            }
            sk++; /* columns index of Mes */
        }
        dg_cell_free(face_cell);
    }

    double vt[Np*Np]; //transpose Vandermonde matrix
    sk=0;
    for(m=0;m<Np;m++){
        for(n=0;n<Np;n++){
            vt[sk++] = dg_cell_V(cell)[n][m];
        }
    }
    double invM[Np*Np];
    matrix_multiply(Np, Np, Np, dg_cell_V(cell)[0], vt, invM);
    /* LIFT = M^{-1}*Mes */
    matrix_multiply(Np, Np, Nfptotal, invM, Mes, LIFT[0]);

    vector_double_free(Mes);
    return LIFT;
}


/**
 * @brief
 * @param cell
 * @param Nfaces
 * @param type
 * @return
 */
dg_cell_face *dg_cell_face_create(dg_cell *cell, const dg_cell_creator *creator){

    dg_cell_face *face = (dg_cell_face *) calloc(1, sizeof(dg_cell_face));
    const int N = dg_cell_N(cell);
    const int Nfaces = dg_cell_Nfaces(cell);
    /* basic info */
    dg_cell_face_info(cell, face);
    /* allocation */
    int *Nfp = face->Nfp;
    double **w = array_double_create(Nfaces, Nfp);
    int **Fmask = array_int_create(Nfaces, Nfp);
    /* assignment of Fmask */
    creator->face_Fmask(cell, Fmask);

//    /* print Fmask */
//    int m,t;
//    for(m=0;m<Nfaces;m++){
//        for(t=0;t<Nfp[m];t++)
//            printf("Fmask[%d][%d]=%d\n", m, t, Fmask[m][t]);
//    }
    /* assignment of ws and Mes */
    dg_cell_type *face_type = dg_cell_facetype(cell);
    int n,f;
    for(f=0;f<Nfaces;f++){
        dg_cell *face_cell = dg_cell_creat(N, face_type[f]);
        int Nface_node = Nfp[f];
        /* assignment of ws */
        for(n=0;n<Nface_node;n++){w[f][n] = dg_cell_w(face_cell)[n];}
        dg_cell_free(face_cell);
    }
    const int Nfptotal = face->Nfptotal;
    face->LIFT = dg_cell_face_LIFT(cell, Nfptotal, Fmask);

    /* assignment */
    face->Fmask = Fmask;
    face->ws = w;
    return face;
}