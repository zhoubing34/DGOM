#include "Convection2d/Convection2d.h"
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
    printf("\n        Deg = %d \n", p_N);
#ifdef Tri
    printf("\n       Triangle Element \n");
#else
    printf("\n     Quadrilateral Element \n");
#endif
    printf("        Ele = %d \n", Ne);
    printf("--------------------------------\n");


#ifdef Tri
    mesh = ReadTriMesh();
    /* find element connections */
    FacePairTri(mesh);

    /* setup mesh */
    SetupTri(mesh);

    // set node connection
    BuildTriMaps(mesh);

    // init mesh coeff
    dt = InitTriMeshInfo(mesh, p_Nfields);
    dt = .5*dt/((p_N+1)*(p_N+1)); // CFL

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


    /* get output */
    Write2TestFile(mesh, FinalTime);
}

void Write2TestFile(Mesh * mesh, double time){

    char filename[10];
    FILE *fig;
    int n, m, std = 0;

    sprintf(filename,"%f",time);
    fig = fopen(filename, "w");

    fprintf(fig, "%s", "\n x=[ \n");
    for (m = 0; m < mesh->K; m ++){
        for (n = 0; n < p_Np; ++n) {
            fprintf(fig, "%f\t", mesh->x[m][n]);
        }
        fprintf(fig, "\n");
    }
    fprintf(fig, "];\n");

    fprintf(fig, "%s", "\n y=[ \n");
    for (m = 0; m < mesh->K; m ++){
        for (n = 0; n < p_Np; ++n) {
            fprintf(fig, "%f\t", mesh->y[m][n]);
        }
        fprintf(fig, "\n");
    }
    fprintf(fig, "];\n");

    fprintf(fig, "%s", "\n var=[ \n");
    for (m = 0; m < mesh->K; m ++){
        for (n = 0; n < p_Np; ++n) {
            fprintf(fig, "%f\t", mesh->f_Q[std++]);
        }
        fprintf(fig, "\n");
    }
    fprintf(fig, "];\n");

    fclose(fig);
}