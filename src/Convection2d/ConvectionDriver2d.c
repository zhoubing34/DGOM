#include <Convection2d/Convection2d.h>
#include <mpi.h>

/*
 * 2d Convection problem
 */

main(int argc, char **argv){

    /* read mesh */
    Mesh * mesh;
    Ncfile * outfile;
    int procid, nprocs, maxNv;
    double dt, FinalTime = 2.4;
    char filename[16] = "Convection2d";

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
    printf("\n    Tri Ele = %d \n", Ne);
#else
    printf("\n   Quad Ele = %d \n", Ne);
#endif
    printf("\n");
    printf("--------------------------------\n");


#ifdef Tri
    mesh = ReadTriMesh();
    /* find element connections */
    FacePairTri(mesh);

    /* setup mesh */
    SetupTriCoeff(mesh);

    // set node connection
    BuildTriMaps(mesh);

    // init mesh coeff
    dt = InitTriMeshInfo(mesh, p_Nfields);
    dt = .5*dt/((p_N+1)*(p_N+1)); // CFL

#else // TODO quad shape
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

    /* setup output file */
    outfile = SetupOutput(mesh, filename);

    /* solve */
    ConvectionRun2d(mesh, outfile ,FinalTime, dt);

    /* TODO output */
//    Write2TestFile(mesh, FinalTime);

    /* finish */
    ConvectionFinish(mesh, outfile);

}

void ConvectionFinish(Mesh * mesh, Ncfile * outfile){

    int ret;

    ret = ncmpi_close(outfile->ncfile);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);


    MPI_Finalize();
}