#include "GenUniformMesh.h"

/**
 * @brief
 * Generation of uniform triangle mesh.
 *
 * @details
 * Generate uniform triangle mesh with specific coordinate range and number of elements.
 * The flag 'type' defines how the square divides into two triangles.
 * The function returns a pointer to `UnstructMesh` object:
 * The number of elements and vertexes - `ne` and `ne`;
 * Variables `EToV`, `vx` and `vy`;
 * All is allocated inside of the function.
 *
 * @param[in]   Mx      number of elements on x coordinate
 * @param[in]   My      number of elements on x coordinate
 * @param[in]   xmin    minimum value of x
 * @param[in]   xmax    maxmum value of x
 * @param[in]   ymin    minimum value of y
 * @param[in]   ymax    maxmum value of y
 * @param[in]   type    flag for triangle dividing, 0 for '\' and 1 for '/'
 *
 * @return      triGrid unstructed uniform mesh of triangle elements
 * @note
 * Generate an uniform triangle grid, the user should call the `DestroyUnstructMesh`
 * function to deallocate the grid after using it.
 */
UnstructMesh* GenUniformTriMesh(int Mx, int My,
                       double xmin, double xmax,
                       double ymin, double ymax, int type){

    /* parameter */
    int Ne = Mx*My*2;
    int Nv = (Mx+1)*(My+1);
    int Dim= 2;

    /* initialize */
    UnstructMesh* grid = CreateUnstructMesh(Dim, Ne, Nv, TRIANGLE);
    grid->name = "uniform triangle grid";

    /* allocation */
    int **EToV = grid->EToV;
    double *VX = grid->vx;
    double *VY = grid->vy;

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

    return grid;
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
 * @param[in]   Mx      number of elements on x coordinate
 * @param[in]   My      number of elements on x coordinate
 * @param[in]   xmin    minimum value of x
 * @param[in]   xmax    maxmum value of x
 * @param[in]   ymin    minimum value of y
 * @param[in]   ymax    maxmum value of y
 *
 * @return      grid    unstructed uniform mesh of quadrilateral elements
 *
 * @note
 * Generate an uniform quadrilateral grid, the user should call the `DestroyUnstructMesh`
 * function to deallocate the grid after using it.
 */
UnstructMesh* GenUniformQuadMesh(int Mx, int My,
                       double xmin, double xmax,
                       double ymin, double ymax){
    /* parameter */
    int Ne = Mx*My;
    int Nv = (Mx+1)*(My+1);
    int Dim=2;

    UnstructMesh* grid = CreateUnstructMesh(Dim, Ne, Nv, QUADRIL);
    grid->name = "uniform quadrilateral grid";

    /* allocation */
    int **EToV = grid->EToV;
    double *VX = grid->vx;
    double *VY = grid->vy;

    /* vertex coordinate */
    int index; /* vertex index */
    int dim1, dim2;
    for (dim1 = 0; dim1<My+1; dim1++){
        for (dim2 = 0; dim2<Mx+1; dim2++){
            /* the order of vertex is from left to right, and from bottom to the top */
            index = dim1*(Mx+1) + dim2;
            VX[index] = (double)dim2/Mx * (xmax - xmin) + xmin;
            VY[index] = (double)dim1/My * (ymax - ymin) + ymin;
        }
    }

    /* EToV connection */
    int ie, topInd[2], botInd[2];
    for (dim1 = 0; dim1 < My; dim1++){
        for (ie = 0; ie < Mx; ie++){
            /* Assignment of vertex index, which start from 0 */
            botInd[0] = dim1*(Mx + 1)+ie;
            botInd[1] = dim1*(Mx + 1)+ie + 1;
            topInd[0] = botInd[0] + Mx + 1;
            topInd[1] = botInd[1] + Mx + 1;

            /* Assignment of EToV,
             * each element start from the left lower vertex */
            EToV[dim1 * Mx + ie][0] = botInd[0];
            EToV[dim1 * Mx + ie][1] = botInd[1];
            EToV[dim1 * Mx + ie][2] = topInd[1];
            EToV[dim1 * Mx + ie][3] = topInd[0];
        }
    }

    return grid;
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
 * @param[in]   Mx      number of elements on x coordinate
 * @param[in]   My      number of elements on x coordinate
 * @param[in]   xmin    minimum value of x
 * @param[in]   xmax    maxmum value of x
 * @param[in]   ymin    minimum value of y
 * @param[in]   ymax    maxmum value of y
 * @param[in]   type    flag for triangle dividing, 0 for '\' and 1 for '/'
 * @param[in]   procid  index of local process
 * @param[in]   nprocs  number of processes
 *
 * @return  parTriGrid  unstructed grid
 *
 * @note
 * Generate a parallel uniform triangle grid, the user should call the `DestroyUnstructMesh`
 * function to deallocate the grid after using it.
 */
UnstructMesh* GenParallelUniformTriMesh(int Mx, int My,
                               double xmin, double xmax,
                               double ymin, double ymax, int type,
                               int procid, int nprocs){

    UnstructMesh *globalGrid = GenUniformTriMesh(Mx, My, xmin, xmax, ymin, ymax, type);
    int K       = globalGrid->ne;
    int **EToV  = globalGrid->EToV;

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

    UnstructMesh *parTriGrid = CreateUnstructMesh(
            globalGrid->dim, Klocal, globalGrid->nv, TRIANGLE);
    parTriGrid->name = "parallel uniform triangle grid";

    /* assignment of vertex coordinate */
    int i;
    for(i=0; i<parTriGrid->nv; i++){
        parTriGrid->vx[i] = globalGrid->vx[i];
        parTriGrid->vy[i] = globalGrid->vy[i];
    }

    /* EToV of local mesh */
    int **newEToV = parTriGrid->EToV;

    /* Assignment of local mesh information */
    int sk=0, n;
    for (n=0; n<K; n++){
        if(n>=Kstart && n<Kstart+Klocal) {
            newEToV[sk][0] = EToV[n][0];
            newEToV[sk][1] = EToV[n][1];
            newEToV[sk][2] = EToV[n][2];
            ++sk;
        }
    }

    free(Kprocs);
    DestroyUnstructMesh(globalGrid);

    return parTriGrid;
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
 * @return  parQuadGrid  unstructed grid
 *
 * @note
 * Generate a parallel uniform quadrilateral grid, the user should call the `DestroyUnstructMesh`
 * function to deallocate the grid after using it.
 */
UnstructMesh* GenParallelUniformQuadMesh(int Mx, int My,
                               double xmin, double xmax,
                               double ymin, double ymax,
                               int procid, int nprocs){

    UnstructMesh *globalGrid = GenUniformQuadMesh(Mx, My, xmin, xmax, ymin, ymax);
    int K       = globalGrid->ne;
    int **EToV  = globalGrid->EToV;

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

    UnstructMesh *parQuadGrid = CreateUnstructMesh(
            globalGrid->dim, Klocal, globalGrid->nv, QUADRIL);
    parQuadGrid->name = "parallel uniform quadrilateral grid";

    /* EToV of local mesh */
    int **newEToV = globalGrid->EToV;

    /* assignment of vertex coordinate */
    int i;
    for(i=0; i<parQuadGrid->nv; i++){
        parQuadGrid->vx[i] = globalGrid->vx[i];
        parQuadGrid->vy[i] = globalGrid->vy[i];
    }

    /* Assignment of local mesh information */
    int sk = 0, n;
    for (n = 0; n < K; n++){
        if(n>=Kstart && n<Kstart+Klocal) {
            newEToV[sk][0] = EToV[n][0];
            newEToV[sk][1] = EToV[n][1];
            newEToV[sk][2] = EToV[n][2];
            newEToV[sk][3] = EToV[n][3];
            ++sk;
        }
    }

    DestroyUnstructMesh(globalGrid);
    free(Kprocs);

    return parQuadGrid;
}