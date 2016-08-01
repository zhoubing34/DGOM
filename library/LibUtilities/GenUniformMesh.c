#include "LibUtilities.h"

/**
 * @brief
 * Generation of uniform triangle mesh.
 *
 * @details
 * Generate uniform triangle mesh with specific coordinate range and number of elements.
 * The flag 'type' defines how the square divides into two triangles.
 * The function returns the number of elements and vertexes - `Ne` and `Nv`.
 * Variables `EToV`, `VX` and `VY` is allocated inside of the function, therefore the pointer is
 * passed as parameters.
 *
 * @param[in] Mx    number of elements on x coordinate
 * @param[in] My    number of elements on x coordinate
 * @param[in] xmin  minimum value of x
 * @param[in] xmax  maxmum value of x
 * @param[in] ymin  minimum value of y
 * @param[in] ymax  maxmum value of y
 * @param[in] type  flag for triangle dividing, 0 for '\' and 1 for '/'
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * newEToV  | int      | vertex index inside of each elements
 * newVX    | double   | coordinate of vertex
 * newVY    | double   | coordinate of vertex
 *
 */
void GenUniformTriMesh(int Mx, int My,
                       double xmin, double xmax,
                       double ymin, double ymax, int type,
                       int *newNe, int *newNv,
                       int ***newEToV, double **newVX, double **newVY){

    /* parameter */
    int Ne = Mx * My *2;
    int Nv = (Mx + 1)*(My + 1);

    *newNe = Ne; *newNv = Nv;

    /* allocation */
    int **EToV = BuildIntMatrix(Ne, 3);
    double *VX = BuildVector(Nv);
    double *VY = BuildVector(Nv);

    /* vertex coordinate */
    int index; /* vertex index */
    int ir, ic;
    for (ir = 0; ir<My+1; ir++){
        for (ic = 0; ic<Mx+1; ic++){
            index = ir*(Mx+1) + ic; /* the order of vertex is from left to right, and from bottom to the top */
            VX[index] = (double)ic/Mx * (xmax - xmin) + xmin;
            VY[index] = (double)ir/My * (ymax - ymin) + ymin;
        }
    }

    /* EToV connection */
    int ie;
    int topNode[2], botNode[2];
    for (ir = 0; ir < My; ir++){
        for (ie = 0; ie < Mx; ie++){
            /* Assignment of vertex index, which start from 0 */
            botNode[0] = ir*(Mx + 1)+ie;        /* bottom left  */
            botNode[1] = ir*(Mx + 1)+ie + 1;    /* bottom right */
            topNode[0] = botNode[0] + (Mx + 1); /* top left  */
            topNode[1] = botNode[1] + (Mx + 1); /* top right */

            /* Assignment of EToV */
            if ( type==1 ){ /* for '/' division */
                EToV[ir * Mx * 2 + ie][0] = botNode[0];
                EToV[ir * Mx * 2 + ie][1] = topNode[1];
                EToV[ir * Mx * 2 + ie][2] = topNode[0];

                EToV[ir * Mx * 2 + Mx + ie][0] = botNode[0];
                EToV[ir * Mx * 2 + Mx + ie][1] = botNode[1];
                EToV[ir * Mx * 2 + Mx + ie][2] = topNode[1];
            }else if( type==0 ){ /* for '\' division */
                EToV[ir * Mx * 2 + ie][0] = botNode[1];
                EToV[ir * Mx * 2 + ie][1] = topNode[1];
                EToV[ir * Mx * 2 + ie][2] = topNode[0];

                EToV[ir * Mx * 2 + Mx + ie][0] = botNode[0];
                EToV[ir * Mx * 2 + Mx + ie][1] = botNode[1];
                EToV[ir * Mx * 2 + Mx + ie][2] = topNode[0];
            }
        }
    }

    /* assignment */
    *newEToV = EToV;
    *newVX   = VX;
    *newVY   = VY;
}

/**
 * @brief
 * Generation of uniform quadrilateral mesh.
 *
 * @details
 * Generate uniform quadrilateral mesh with specific coordinate range and number of elements.
 * The function returns the number of elements and vertexes - `Ne` and `Nv`.
 * Variables `EToV`, `VX` and `VY` is allocated inside of the function, therefore the pointer is
 * passed as parameters.
 *
 * @param[in] Mx    number of elements on x coordinate
 * @param[in] My    number of elements on x coordinate
 * @param[in] xmin  minimum value of x
 * @param[in] xmax  maxmum value of x
 * @param[in] ymin  minimum value of y
 * @param[in] ymax  maxmum value of y
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * newEToV  | int      | vertex index inside of each elements
 * newVX    | double   | coordinate of vertex
 * newVY    | double   | coordinate of vertex
 *
 */
void GenUniformQuadMesh(int Mx, int My,
                       double xmin, double xmax,
                       double ymin, double ymax,
                       int *newNe, int *newNv,
                       int ***newEToV, double **newVX, double **newVY){
    /* parameter */
    int Ne = Mx * My;
    int Nv = (Mx + 1)*(My + 1);

    *newNe = Ne; *newNv = Nv;

    /* allocation */
    int **EToV = BuildIntMatrix(Ne, 4);
    double *VX = BuildVector(Nv);
    double *VY = BuildVector(Nv);

    /* vertex coordinate */
    int index; /* vertex index */
    int ir, ic;
    for (ir = 0; ir<My+1; ir++){
        for (ic = 0; ic<Mx+1; ic++){
            index = ir*(Mx+1) + ic; /* the order of vertex is from left to right, and from bottom to the top */
            VX[index] = (double)ic/Mx * (xmax - xmin) + xmin;
            VY[index] = (double)ir/My * (ymax - ymin) + ymin;
        }
    }

    /* EToV connection */
    int ie;
    int topNode[2], botNode[2];
    for (ir = 0; ir < My; ir++){
        for (ie = 0; ie < Mx; ie++){
            /* Assignment of vertex index, which start from 0 */
            botNode[0] = ir*(Mx + 1)+ie;
            botNode[1] = ir*(Mx + 1)+ie + 1;
            topNode[0] = botNode[0] + Mx + 1;
            topNode[1] = botNode[1] + Mx + 1;

            /* Assignment of EToV,
             * each element start from the left lower vertex */
            EToV[ir * Mx + ie][0] = botNode[0];
            EToV[ir * Mx + ie][1] = botNode[1];
            EToV[ir * Mx + ie][2] = topNode[1];
            EToV[ir * Mx + ie][3] = topNode[0];
        }
    }

    /* assignment */
    *newEToV = EToV;
    *newVX   = VX;
    *newVY   = VY;

}

/**
 * @brief
 * Generate uniform triangle mesh and separate it equally on multiply processes.
 *
 * @details
 * Generate uniform triangle mesh with specific coordinate range and number of elements on each edge.
 * The flag `type` defines how the square divides into two triangles.
 * The variable `parEToV` only contains the local elements distributed on process `procid`.
 * The function returns the number of elements and vertexes - `parK` and `newNv`.
 * The vertex is not separated and all the coordinates is stored in `VX` and `VY`.
 * Variables `parEToV`, `VX` and `VY` is allocated inside of the function, therefore the pointer is
 * passed as parameters.
 *
 * @param[in] Mx     number of elements on x coordinate
 * @param[in] My     number of elements on x coordinate
 * @param[in] xmin   minimum value of x
 * @param[in] xmax   maxmum value of x
 * @param[in] ymin   minimum value of y
 * @param[in] ymax   maxmum value of y
 * @param[in] type   flag for triangle dividing, 0 for '\' and 1 for '/'
 * @param[in] procid index of process
 * @param[in] nprocs number of processes
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * parEToV  | int      | vertex index inside of local elements
 * VX    | double   | coordinate of vertex
 * VY    | double   | coordinate of vertex
 *
 */
void GenParallelUniformTriMesh(int Mx, int My,
                               double xmin, double xmax,
                               double ymin, double ymax, int type,
                               int procid, int nprocs,
                               int *parK, int *newNv,
                               int ***parEToV, double **VX, double **VY){

    int **EToV, Ne, Nv;
    GenUniformTriMesh(Mx, My, xmin, xmax, ymin, ymax, type, &Ne, &Nv, &EToV, VX, VY);
    int K       = Ne;
    int Nvertex = 3;

    int *Kprocs = (int *) calloc(nprocs, sizeof(int));
    int  Klocal = (int)( (double)K/ (double)nprocs );

    /* number of elements in each process */
    int p;
    for(p = 0; p < nprocs - 1; ++p){
        Kprocs[p] = Klocal;
    }
    /* the rest elements is distributed to the last process: p = nprocs - 1  */
    Kprocs[p] = K + Klocal - nprocs*Klocal;

    /* number of element in local process */
    Klocal = Kprocs[procid];

    /* start index of element */
    int Kstart = 0;
    for(p=0;p<procid;++p)
        Kstart += Kprocs[p];


    /* EToV of local mesh */
    int **newEToV = BuildIntMatrix(Klocal, Nvertex);
    /* Assignment of local mesh information */
    int sk = 0, ie;
    for (ie = 0; ie < K; ie++){
        if(ie>=Kstart && ie<Kstart+Klocal) {
            newEToV[sk][0] = EToV[ie][0];
            newEToV[sk][1] = EToV[ie][1];
            newEToV[sk][2] = EToV[ie][2];
            ++sk;
        }
    }

    DestroyIntMatrix(EToV);

    *parEToV = newEToV;
    *parK    = Klocal;
    *newNv   = Nv;
}

/**
 * @brief
 * Generation of uniform quadrilateral mesh.
 *
 * @details
 * Generate uniform quadrilateral mesh with specific coordinate range and number of elements.
 * The variable `parEToV` only contains the local elements distributed on process `procid`.
 * The function returns the number of elements and vertexes - `parK` and `newNv`.
 * The vertex is not separated and all the coordinates is stored in `VX` and `VY`.
 * Variables `parEToV`, `VX` and `VY` is allocated inside of the function, therefore the pointer is
 * passed as parameters.
 *
 * @param[in] Mx    number of elements on x coordinate
 * @param[in] My    number of elements on x coordinate
 * @param[in] xmin  minimum value of x
 * @param[in] xmax  maxmum value of x
 * @param[in] ymin  minimum value of y
 * @param[in] ymax  maxmum value of y
 * @param[in] procid index of process
 * @param[in] nprocs number of processes
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * parEToV  | int      | vertex index inside of local elements
 * VX    | double   | coordinate of vertex
 * VY    | double   | coordinate of vertex
 *
 */
void GenParallelUniformQuadMesh(int Mx, int My,
                               double xmin, double xmax,
                               double ymin, double ymax,
                               int procid, int nprocs,
                               int *parK, int *newNv,
                               int ***parEToV, double **VX, double **VY){

    int **EToV, Ne, Nv;
    GenUniformQuadMesh(Mx, My, xmin, xmax, ymin, ymax, &Ne, &Nv, &EToV, VX, VY);
    int K       = Ne;
    int Nvertex = 4;

    int *Kprocs = (int *) calloc(nprocs, sizeof(int));
    int  Klocal = (int)( (double)K/ (double)nprocs );

    /* number of elements in each process */
    int p;
    for(p = 0; p < nprocs - 1; ++p){
        Kprocs[p] = Klocal;
    }
    /* the rest elements is distributed to the last process: p = nprocs - 1  */
    Kprocs[p] = K + Klocal - nprocs*Klocal;

    /* number of element in local process */
    Klocal = Kprocs[procid];

    /* start index of element */
    int Kstart = 0;
    for(p=0;p<procid;++p)
        Kstart += Kprocs[p];


    /* EToV of local mesh */
    int **newEToV = BuildIntMatrix(Klocal, Nvertex);
    /* Assignment of local mesh information */
    int sk = 0, ie;
    for (ie = 0; ie < K; ie++){
        if(ie>=Kstart && ie<Kstart+Klocal) {
            newEToV[sk][0] = EToV[ie][0];
            newEToV[sk][1] = EToV[ie][1];
            newEToV[sk][2] = EToV[ie][2];
            newEToV[sk][3] = EToV[ie][3];
            ++sk;
        }
    }

    DestroyIntMatrix(EToV);

    *parEToV = newEToV;
    *parK    = Klocal;
    *newNv   = Nv;
}