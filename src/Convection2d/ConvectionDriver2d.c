/** \file   ConvectionDriver2d.c
 *  \brief  2d scalar convection problem
 */

#include <Convection2d/Convection2d.h>
#include <mpi.h>

/**
 * @brief
 * main function for 2d convection problem
 *
 * @details
 * 2d scalar convection problem
 * \f[ \frac{\partial C}{\partial t} +
 * \frac{\partial uC}{\partial x} + \frac{\partial vC}{\partial y} = 0 \f]
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 * @note
 * The Mass, Derivative and LIFT matrix for each order is include from the header files.
 * In the future, these matrix will be generate from library.
 *
 * @todo
 * 1. quadrilateral shapes
 * 2. slope limiter for high order shapes
 */
main(int argc, char **argv){

    /* read mesh */
    Mesh * mesh;
    Ncfile * outfile;
    int procid, nprocs;
    double dt, FinalTime = 2.4;
    char filename[16] = "Convection2d";

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if(!procid) {
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
    }

#if defined TRI
    mesh = ReadTriMesh();

#if 0 // check mesh
    PrintMeshTri(mesh);
#endif

    /* redistribute element */
    LoadBalanceTri(mesh);

    /* find element connections */
    FacePairTri(mesh);

    /* setup mesh */
    SetupTriCoeff(mesh);

    /* set node connection */
    BuildTriMaps(mesh);

    /* init mesh coeff */
    dt = InitTriMeshInfo(mesh, p_Nfields);
    dt = .5*dt/((p_N+1)*(p_N+1)); // CFL

#elif defined QUAD // TODO quad shape
    mesh = ReadQuadMesh();

    /* find element connections */
    FacePairQuad(mesh);
#endif

    // check mesh connection
#if 0
    PrintMeshConnectionTri(mesh);
#endif

    /* initial conditions */
    InitData(mesh);

    /* setup output file */
    outfile = SetupOutput(mesh, filename);

    /* solve */
    ConvectionRun2d(mesh, outfile ,FinalTime, dt);

    /* finish */
    ConvectionFinish(mesh, outfile);

}


void ConvectionFinish(Mesh * mesh, Ncfile * outfile){

    int ret;
    ret = ncmpi_close(outfile->ncfile);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    MPI_Finalize();
}