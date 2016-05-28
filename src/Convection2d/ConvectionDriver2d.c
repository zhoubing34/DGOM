#include "Convection2d/ConvectionDriver.h"
#include <mpi.h>

/*
 * 2d scalar convection problem
 */

main(int argc, char **argv){

    /* read mesh */
    Mesh * mesh;
    int procid, nprocs, maxNv;
    double dt, FinalTime = 2.4;

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    /* assign gpu */
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


    printf("--------------------------------\n");
    printf("          Convection2d\n");
    printf("--------------------------------\n");
    printf("\n    2d Convection Test Case\n");
    printf("\n        Node = %d \n", p_Np);
    printf("\n        Ele  = %d \n", Ne);
    printf("--------------------------------\n");


#ifdef Tri
    mesh = ReadTriMesh();
    /* find element connections */
    FacePairTri(mesh);

    /* setup mesh */
    SetupTri(mesh);

    // set node connection
    BuildTriMaps(mesh);

    // init mesh coefficent
    dt = InitTriMeshInfo(mesh, p_Nfields);


#else
    mesh = ReadQuadMesh();

    /* find element connections */
    FacePairQuad(mesh);

#endif

#if 0 // check mesh
    PrintMeshTri(mesh);
#endif

#if 0 // check mesh connection
    PrintMeshConnectionTri(mesh);
#endif


    /* initial conditions */
    InitData(mesh);

    /* solve */
    ConvectionRun2d(mesh, FinalTime, dt);

    /* get max and min */

}