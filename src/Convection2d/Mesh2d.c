#include "Convection2d/Convection2d.h"
#include <mpi.h>

/*
 * Set local mesh grid
 * */
Mesh* ReadTriMesh(){
    // generate uniform triangle mesh
    // ne - # of element on each edge

    Mesh *mesh = (Mesh*) calloc(1, sizeof(Mesh));
    int upIndex[2], downIndex[2]; // up and down layer vertex index
    int irow, icol, p, ie, K, Klocal, Kstart;
    int * Kprocs;

    /* decide on parition */
    int procid, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if(!procid) printf("Root: Entering ReadTriMesh\n");

    mesh->procid = procid;
    mesh->nprocs = nprocs;

    mesh->Nverts = 3; /* assume triangles */
    mesh->Nfaces = 3; /* assume triangles */

    Kprocs = (int *) calloc(nprocs, sizeof(int));
    K = Ne* Ne *2; /* total # of elements */
    Klocal = (int)( (double)K/ (double)nprocs );

    /* # of elements in each process */
    for(p = 0; p < nprocs - 1; ++p){
        Kprocs[p] = Klocal;
    }
    Kprocs[p] = K + Klocal - nprocs*Klocal;

    Klocal = Kprocs[procid]; // local process element number

    Kstart= 0; // local element start index
    for(p=0;p<procid;++p)
        Kstart += Kprocs[p];

    mesh->Nv = (Ne + 1)*(Ne + 1); // # of vertex

    // allocate vertex memory
    double *VX = BuildVector(mesh->Nv);
    double *VY = BuildVector(mesh->Nv);

    // allocate element memory
    int **newEToV = BuildIntMatrix(K, mesh->Nverts);
    mesh->EToV = BuildIntMatrix(Klocal, mesh->Nverts);

    for (irow = 0; irow < Ne; irow++){
        for (ie = 0; ie < Ne; ie++){
            // vertex index, start from 0
            upIndex[0] = irow*(Ne + 1)+ie;
            upIndex[1] = irow*(Ne + 1)+ ie + 1;
            downIndex[0] = upIndex[0] + Ne + 1;
            downIndex[1] = upIndex[1] + Ne + 1;

            // set EToV
#define UP_RIGHT
#ifdef UP_RIGHT
            newEToV[irow * Ne * 2 + ie][0] = downIndex[0];
            newEToV[irow * Ne * 2 + ie][1] = upIndex[1];
            newEToV[irow * Ne * 2 + ie][2] = upIndex[0];

            newEToV[irow * Ne * 2 + Ne + ie][0] = downIndex[0];
            newEToV[irow * Ne * 2 + Ne + ie][1] = downIndex[1];
            newEToV[irow * Ne * 2 + Ne + ie][2] = upIndex[1];
#else // up-left
            mesh->EToV[irow * Ne * 2 + ie][0] = downIndex[1];
            mesh->EToV[irow * Ne * 2 + ie][1] = upIndex[1];
            mesh->EToV[irow * Ne * 2 + ie][2] = upIndex[0];

            mesh->EToV[irow * Ne * 2 + Ne + ie][0] = downIndex[0];
            mesh->EToV[irow * Ne * 2 + Ne + ie][1] = downIndex[1];
            mesh->EToV[irow * Ne * 2 + Ne + ie][2] = upIndex[0];
#endif
        }
    }

    int index;

    for (irow = 0; irow<Ne+1; irow++){
        for (icol = 0; icol<Ne+1; icol++){
            index = irow*(Ne+1) + icol;
            VX[index] = (double)icol/Ne *2.0 -1;
            VY[index] = (double)irow/Ne * (-2.0) + 1;
        }
    }

    mesh->GX = BuildMatrix(Klocal, mesh->Nverts);
    mesh->GY = BuildMatrix(Klocal, mesh->Nverts);

    int sk = 0;
    for (ie = 0; ie < K; ie++){
        if(ie>=Kstart && ie<Kstart+Klocal) {
            mesh->EToV[sk][0] = newEToV[ie][0];
            mesh->EToV[sk][1] = newEToV[ie][1];
            mesh->EToV[sk][2] = newEToV[ie][2];

            mesh->GX[sk][0] = VX[mesh->EToV[sk][0]];
            mesh->GX[sk][1] = VX[mesh->EToV[sk][1]];
            mesh->GX[sk][2] = VX[mesh->EToV[sk][2]];

            mesh->GY[sk][0] = VY[mesh->EToV[sk][0]];
            mesh->GY[sk][1] = VY[mesh->EToV[sk][1]];
            mesh->GY[sk][2] = VY[mesh->EToV[sk][2]];
            ++sk;

        }
    }

    mesh->K = Klocal;

    DestroyVector(VX);
    DestroyVector(VY);
    DestroyIntMatrix(newEToV);

    if(!procid) printf("Root: Leaving ReadTriMesh\n");
    return mesh;
}

Mesh* ReadQuadMesh(){
    // generate uniform quad mesh

    Mesh * mesh;
    int i;

    mesh->Nverts = 4; /* assume quadrial */
    mesh->Nfaces = 4; /* assume quadrial */

    mesh->K = Ne* Ne *2; // # of elements
    mesh->Nv = (Ne + 1)*(Ne + 1); // # of vertex


    double *VX = BuildVector(mesh->Nv);
    double *VY = BuildVector(mesh->Nv);

    for (i = 0; i < Ne+1; i++){

    }

    DestroyVector(VX);
    DestroyVector(VY);

    return mesh;

}

void GeometricFactorsTri(Mesh *mesh, int k,
                        double *drdx, double *dsdx, double *drdy, double *dsdy,
                        double *J){

    double x1 = mesh->GX[k][0], y1 =  mesh->GY[k][0];
    double x2 = mesh->GX[k][1], y2 =  mesh->GY[k][1];
    double x3 = mesh->GX[k][2], y3 =  mesh->GY[k][2];

    /* compute geometric factors of the following afine map */
    /* x = 0.5*( -(r+s)*x1 + (1+r)*x2 + (1+s)*x3 ) */
    /* y = 0.5*( -(r+s)*y1 + (1+r)*y2 + (1+s)*y3 ) */

    double dxdr = (x2-x1)/2,  dxds = (x3-x1)/2;
    double dydr = (y2-y1)/2,  dyds = (y3-y1)/2;

    /* Jacobian of coordinate mapping */
    *J = -dxds*dydr + dxdr*dyds;

    if(*J<0)
        printf("warning: J = %lg\n", *J);

    /* inverted Jacobian matrix for coordinate mapping */
    *drdx =  dyds/(*J);
    *dsdx = -dydr/(*J);
    *drdy = -dxds/(*J);
    *dsdy =  dxdr/(*J);
}

void NormalsTri(Mesh *mesh, int k, double *nx, double *ny, double *sJ){

    int f;

    double x1 = mesh->GX[k][0], y1 = mesh->GY[k][0];
    double x2 = mesh->GX[k][1], y2 = mesh->GY[k][1];
    double x3 = mesh->GX[k][2], y3 = mesh->GY[k][2];

    nx[0] =  (y2-y1);  ny[0] = -(x2-x1);
    nx[1] =  (y3-y2);  ny[1] = -(x3-x2);
    nx[2] =  (y1-y3);  ny[2] = -(x1-x3);

    for(f=0;f<3;++f){
        sJ[f] = sqrt(nx[f]*nx[f]+ny[f]*ny[f]);
        nx[f] /= sJ[f];
        ny[f] /= sJ[f];
        sJ[f] /= 2.;
    }
}

#define DSET_NAME_LEN 1024

void PrintMeshTri ( Mesh *mesh ){
    int n, m, rank, nprocs, ret;
    char filename[DSET_NAME_LEN];

    rank = mesh->procid;
    nprocs = mesh->nprocs;

    ret = snprintf(filename, DSET_NAME_LEN, "%d-%d.txt", rank, nprocs);
    if (ret >= DSET_NAME_LEN) {
        fprintf(stderr, "name too long \n");
        exit(-1);
    }

    FILE *fig = fopen(filename, "w");

    fprintf(fig, "Mesh data: \n");
    fprintf(fig, "\n K = %d\n", mesh->K);
    fprintf(fig, "\n Nv = %d\n", mesh->Nv);
    fprintf(fig, "\n Nverts = %d\n", mesh->Nverts);
    fprintf(fig, "\n Vertex coordinates X: \n");
    for (n = 0; n < mesh->K; n ++){
        fprintf(fig, "%f, \t %f, \t %f\n", mesh->GX[n][0], mesh->GX[n][1], mesh->GX[n][2]);
    }

    fprintf(fig, "\n Vertex coordinates Y: \n");
    for (n = 0; n < mesh->K; n ++){
        fprintf(fig, "%f, \t %f, \t %f\n", mesh->GY[n][0], mesh->GY[n][1], mesh->GY[n][2]);
    }

    fprintf(fig, "\n Element to vertex connectivity = \n");
    for(n=0;n<mesh->K;++n){
        fprintf(fig, "%d: %d %d %d \n", n,
               mesh->EToV[n][0], mesh->EToV[n][1], mesh->EToV[n][2]);
    }

//    fprintf(fig, "\n Node coordinates: \n");
//    for (m = 0; m<mesh->K; m++){
//        for (n = 0; n<p_Np; n++)
//            fprintf(fig, "[%f, %f] \t", mesh->x[m][n], mesh->y[m][n]);
//        fprintf(fig, "\n");
//    }
    fclose(fig);

}

