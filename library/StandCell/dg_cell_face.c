//
// Created by li12242 on 17/3/14.
//

#include "dg_cell_face.h"

/**
 * @brief deallocate the memory of dg_cell_face type pointer.
 * @param face pointer of dg_cell_face type;
 */
void dg_cell_face_free(dg_cell_face *face){
    vector_int_free(face->Nfp);
    array_double_free(face->ws);
    array_int_free(face->Fmask);
    matrix_double_free(face->LIFT);
    free(face);
    return;
}

/**
 * @brief set the properties of Nfptotal and Nfp in face.
 * @param cell pointer of dg_cell type;
 * @param face pointer of dg_cell_face type;
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
/**
 * @brief Calculate the LIFT matrix.
 * @details
 * The LIFT matrix is to lift the surface integral flux into nodal terms.
 * @param cell pointer of dg_cell type;
 * @param Nfptotal total number of nodes on each faces;
 * @param Fmask nodex index on all faces;
 * @return
 * pointer to a new LIFT matrix.
 */
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
 * @brief calculate the Fmask matrix.
 * @details
 * Fmask matrix contains all the nodal index on faces. This program searches
 * all the nodal points to match all the face nodes.
 * @param cell pointer of dg_cell type;
 * @param Fmask the Fmask matrix;
 */
void dg_face_cell_Fmask(dg_cell *cell, int **Fmask){
    const int N = dg_cell_N(cell);
    const int Np = dg_cell_Np(cell);
    const int Nfaces = dg_cell_Nfaces(cell);
    const dg_cell_type *face_type = dg_cell_facetype(cell);
    int m,n,f;
    for(f=0;f<Nfaces;f++){
        dg_cell *face_cell = dg_cell_creat(N, face_type[f]);
        /* assignment of Fmask */
        int Nface_node = dg_cell_Np(face_cell);
        int Nface_vert = dg_cell_Nv(face_cell);
        double f_vr[Nface_vert], f_vs[Nface_vert], f_vt[Nface_vert];
        for(n=0;n<Nface_vert;n++){
            f_vr[n] = dg_cell_vr(cell)[ dg_cell_FToV(cell)[f][n] ];
            f_vs[n] = dg_cell_vs(cell)[ dg_cell_FToV(cell)[f][n] ];
            f_vt[n] = dg_cell_vt(cell)[ dg_cell_FToV(cell)[f][n] ];
        }
        // face node coordinate
        double f_r[Nface_node], f_s[Nface_node], f_t[Nface_node];
        // get the face node coordinate by map from vertex
        face_cell->proj_vert2node(face_cell, f_vr, f_r);
        face_cell->proj_vert2node(face_cell, f_vs, f_s);
        face_cell->proj_vert2node(face_cell, f_vt, f_t);
        for(n=0;n<Nface_node;n++){
            double x1 = f_r[n];
            double y1 = f_s[n];
            double z1 = f_t[n];
            for(m=0;m<Np;m++){ // loop over all the nodes
                double x2 = dg_cell_r(cell)[m];
                double y2 = dg_cell_s(cell)[m];
                double z2 = dg_cell_t(cell)[m];
                double d12 = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
                if(d12 < EPS) {Fmask[f][n] = m; break;}
            }
        }
    }
    return;
}

/**
 * @brief creating the dg_cell_face structure.
 * @param cell pointer of dg_cell type;
 * @return
 * pointer to a new dg_cell_face structure.
 */
dg_cell_face *dg_cell_face_create(dg_cell *cell){
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
    dg_face_cell_Fmask(cell, Fmask);

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