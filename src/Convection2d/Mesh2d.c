#include "ConvectionDriver2d.h"

/* private function */
MultiReg2d* ReadTriMesh(StdRegions2d *shape, int Ne);
MultiReg2d* ReadQuadMesh(StdRegions2d *shape, int Ne);

#define UP_RIGHT // define triangle shapes

MultiReg2d* ReadMesh(StdRegions2d *shape, int Ne){
    MultiReg2d *mesh;

    if (shape->Nfaces == 3){
        mesh = ReadTriMesh(shape, Ne);
    }else if(shape->Nfaces == 4){
        mesh = ReadQuadMesh(shape, Ne);
    }else{
        fprintf(stderr, "Wrong number of element shapes, Nfaces=%d \n", shape->Nfaces);
        exit(-1);
    }
    return mesh;
}

/**
 * @brief
 * Generate local uniform triangle mesh for each processor
 *
 * @details
 * The whole domain is [-1, 1]x[-1, 1], with each edge partition to *Ne* 
 * elements. The vertex list is arranged from upper left to lower right
 * directions. The whole mesh is generate with newEToV for element to 
 * vertex list. Each process reads part of mesh into mesh->EToV and its 
 * corresponding vertex coordinate mesh->GX and mesh->GY. The uniform 
 * triangle element is separated from the square, and the separation is 
 * controlled by UP_RIGHT macro.
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 * @return
 * return values：
 * name     | type      | description of value
 * -------- |---------- |----------------------
 * mesh     | Mesh*     | mesh object
 *
 */
MultiReg2d* ReadTriMesh(StdRegions2d *shape, int Ne){

    /* mesh object and allocation */
    MultiReg2d *mesh;
    /* index of vertex at upper and lower layer */
    int upperNode[2], lowerNode[2];
    int irow, icol, ie, p, K;
    /** number of elements in local processor */
    int Klocal;
    /** start index of element in all element */
    int Kstart = 0;
    /** number of elements in each processor */
    int *Kprocs;

    /* decide on parition */
    int procid, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    Kprocs = (int *) calloc(nprocs, sizeof(int));
    K = Ne* Ne *2; /* total # of elements */
    Klocal = (int)( (double)K/ (double)nprocs );

    /* number of elements in each process */
    for(p = 0; p < nprocs - 1; ++p){
        Kprocs[p] = Klocal;
    }
    /* the rest elements is distributed to the last process: p = nprocs - 1  */
    Kprocs[p] = K + Klocal - nprocs*Klocal;

    /* number of element in local process */
    Klocal = Kprocs[procid];

    /* start index of element */
    for(p=0;p<procid;++p)
        Kstart += Kprocs[p];

    /* element to vertex list for all elements */
    int **newEToV = BuildIntMatrix(K, shape->Nv);
    /* EToV of local mesh */

    for (irow = 0; irow < Ne; irow++){
        for (ie = 0; ie < Ne; ie++){
            /* Assignment of vertex index, which start from 0 */
            upperNode[0] = irow*(Ne + 1)+ie;
            upperNode[1] = irow*(Ne + 1)+ ie + 1;
            lowerNode[0] = upperNode[0] + Ne + 1;
            lowerNode[1] = upperNode[1] + Ne + 1;

            /* Assignment of EToV */
#ifdef UP_RIGHT
            newEToV[irow * Ne * 2 + ie][0] = lowerNode[0];
            newEToV[irow * Ne * 2 + ie][1] = upperNode[1];
            newEToV[irow * Ne * 2 + ie][2] = upperNode[0];

            newEToV[irow * Ne * 2 + Ne + ie][0] = lowerNode[0];
            newEToV[irow * Ne * 2 + Ne + ie][1] = lowerNode[1];
            newEToV[irow * Ne * 2 + Ne + ie][2] = upperNode[1];
#else // up-left
            newEToV[irow * Ne * 2 + ie][0] = lowerNode[1];
            newEToV[irow * Ne * 2 + ie][1] = upperNode[1];
            newEToV[irow * Ne * 2 + ie][2] = upperNode[0];

            newEToV[irow * Ne * 2 + Ne + ie][0] = lowerNode[0];
            newEToV[irow * Ne * 2 + Ne + ie][1] = lowerNode[1];
            newEToV[irow * Ne * 2 + Ne + ie][2] = upperNode[0];
#endif
        }
    }

    /* # of all vertex */
    int Nverts = (Ne + 1)*(Ne + 1);

    /* allocate coordinate of all vertex */
    double *VX = BuildVector(Nverts);
    double *VY = BuildVector(Nverts);

    /* Assignment of coordinate for vertex */
    int index; /* vertex index */
    for (irow = 0; irow<Ne+1; irow++){
        for (icol = 0; icol<Ne+1; icol++){
            index = irow*(Ne+1) + icol;
            VX[index] = (double)icol/Ne *2.0 -1;
            VY[index] = (double)irow/Ne * (-2.0) + 1;
        }
    }

    int **EToV = BuildIntMatrix(Klocal, shape->Nv);

    /* Assignment of local mesh information */
    int sk = 0;
    for (ie = 0; ie < K; ie++){
        if(ie>=Kstart && ie<Kstart+Klocal) {
            EToV[sk][0] = newEToV[ie][0];
            EToV[sk][1] = newEToV[ie][1];
            EToV[sk][2] = newEToV[ie][2];
            ++sk;

        }
    }

    /* generate triangle mesh */
    mesh = GenMultiReg2d(shape, Klocal, Nverts, EToV, VX, VY);

    DestroyVector(VX);
    DestroyVector(VY);
    DestroyIntMatrix(newEToV);
    DestroyIntMatrix(EToV);

    free(Kprocs);
    return mesh;
}

/**
 * @brief
 * Generate uniform quadrilateral mesh
 *
 * @details
 * The whole domain is [-1, 1]x[-1, 1], with each edge partition to *Ne* 
 * elements. The vertex list is arranged from upper left to lower right
 * directions. The whole mesh is generate with newEToV for element to 
 * vertex list. Each processor reads part of mesh into mesh->EToV and 
 * its corresponding vertex coordinate mesh->GX and mesh->GY.
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 * @return
 * return values：
 * name     | type      | description of value
 * -------- |---------- |----------------------
 * mesh     | Mesh*     | mesh object
 *
 */
MultiReg2d* ReadQuadMesh(StdRegions2d *shape, int Ne){
    /* mesh object and allocation */
    MultiReg2d *mesh;
    /* index of vertex at upper and lower layer */
    int upperNode[2], lowerNode[2];
    int irow, icol, ie, p, K;
    /** number of elements in local processor */
    int Klocal;
    /** start index of element in all element */
    int Kstart = 0;
    /** number of elements in each processor */
    int * Kprocs;

    /* decide on parition */
    int procid, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    Kprocs = (int *) calloc(nprocs, sizeof(int));
    K = Ne* Ne; /* total # of elements */
    Klocal = (int)( (double)K/ (double)nprocs );

    /* number of elements in each process */
    for(p = 0; p < nprocs - 1; ++p){
        Kprocs[p] = Klocal;
    }
    /* the rest elements is distributed to the last process: p = nprocs - 1  */
    Kprocs[p] = K + Klocal - nprocs*Klocal;

    /* number of element in local process */
    Klocal = Kprocs[procid];

    /* start index of element */
    for(p=0;p<procid;++p)
        Kstart += Kprocs[p];

    /* # of all vertex */
    int Nverts = (Ne + 1)*(Ne + 1);

    /* allocate coordinate of all vertex */
    double *VX = BuildVector(Nverts);
    double *VY = BuildVector(Nverts);

    /* element to vertex list for all elements */
    int **newEToV = BuildIntMatrix(K, shape->Nv);
    /* EToV of local mesh */
    int **EToV = BuildIntMatrix(Klocal, shape->Nv);

    for (irow = 0; irow < Ne; irow++){
        for (ie = 0; ie < Ne; ie++){
            /* Assignment of vertex index, which start from 0 */
            upperNode[0] = irow*(Ne + 1)+ie;
            upperNode[1] = irow*(Ne + 1)+ ie + 1;
            lowerNode[0] = upperNode[0] + Ne + 1;
            lowerNode[1] = upperNode[1] + Ne + 1;

            /* Assignment of EToV,
             * each element start from the left lower vertex */
            newEToV[irow * Ne + ie][0] = lowerNode[0];
            newEToV[irow * Ne + ie][1] = lowerNode[1];
            newEToV[irow * Ne + ie][2] = upperNode[1];
            newEToV[irow * Ne + ie][3] = upperNode[0];
        }
    }

    /* Assignment of coordinate for vertex */
    int index; /* vertex index */
    for (irow = 0; irow<Ne+1; irow++){
        for (icol = 0; icol<Ne+1; icol++){
            index = irow*(Ne+1) + icol;
            VX[index] = (double)icol/Ne *2.0 -1;
            VY[index] = (double)irow/Ne * (-2.0) + 1;
        }
    }

    /* Assignment of local mesh information */
    int sk = 0;
    for (ie = 0; ie < K; ie++){
        if(ie>=Kstart && ie<Kstart+Klocal) {
            EToV[sk][0] = newEToV[ie][0];
            EToV[sk][1] = newEToV[ie][1];
            EToV[sk][2] = newEToV[ie][2];
            EToV[sk][3] = newEToV[ie][3];
            ++sk;
        }
    }

    mesh = GenMultiReg2d(shape, Klocal, Nverts, EToV, VX, VY);

    /* deallocate memory */
    DestroyVector(VX);
    DestroyVector(VY);
    DestroyIntMatrix(newEToV);
    DestroyIntMatrix(EToV);

    free(Kprocs);
    return mesh;
}


/**
 * @brief
 * Print the mesh information of each process in files
 *
 * @details
 *
 * @author 
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 */
#define DSET_NAME_LEN 1024
void PrintMesh ( MultiReg2d *mesh ){
    int n, m, rank, nprocs, ret, sk, sf;
    char filename[DSET_NAME_LEN];

    StdRegions2d *shape = mesh->stdcell;

    rank = mesh->procid;
    nprocs = mesh->nprocs;

    ret = snprintf(filename, DSET_NAME_LEN, "%s.%d-%d.txt", "MeshCheck",rank, nprocs);
    if (ret >= DSET_NAME_LEN) {
        fprintf(stderr, "name too long \n");
        exit(-1);
    }

    FILE *fig = fopen(filename, "w");

    fprintf(fig, "Mesh data: \n");
    fprintf(fig, "\n K = %d\n", mesh->K);
    fprintf(fig, "\n Nv = %d\n", shape->Nv);
    fprintf(fig, "\n Total Nverts = %d\n", mesh->Nv);
    fprintf(fig, "\n Vertex coordinates X: \n");
    for (n = 0; n < mesh->K; n ++){
        for( m = 0; m < shape->Nv; m++)
            fprintf(fig, "%f,\t", mesh->GX[n][m]);
        fprintf(fig, "\n");
    }

    fprintf(fig, "\n Vertex coordinates Y: \n");
    for (n = 0; n < mesh->K; n ++){
        for( m = 0; m < shape->Nv; m++)
            fprintf(fig, "%f,\t", mesh->GY[n][m]);
        fprintf(fig, "\n");
    }

#if 1
    fprintf(fig, "\n Node coordinate = \n");
    for (n = 0; n < mesh->K; n ++){
        for( m = 0; m < shape->Np; m++)
            fprintf(fig, "[%f, %f],\t", mesh->x[n][m], mesh->y[n][m]);
        fprintf(fig, "\n");
    }

#endif

#if 1

    sk = 0;
    fprintf(fig, "\n Node vgeo = \n");
    for (n = 0; n < mesh->K; n ++){
        fprintf(fig, "Ele %d:\n", n);
        for( m = 0; m < shape->Np; m++)
            fprintf(fig, "drdx = %f, drdy = %f, dsdx = %f, dsdy = %f\n",
                    mesh->vgeo[sk++], mesh->vgeo[sk++], mesh->vgeo[sk++], mesh->vgeo[sk++]);
        fprintf(fig, "\n");
    }
#endif
    fclose(fig);

}

