/**
 * @file
 * ConvectionDriver2d.c
 *
 * @brief
 * 2D scalar convection problem main function
 */

#include <Convection2d/Convection2d.h>

/**
 * @brief
 * main function for 2d convection problem
 *
 * @details
 * 2d scalar convection problem
 * \f[ \frac{\partial C}{\partial t} + \frac{\partial uC}{\partial x}
 * + \frac{\partial vC}{\partial y} = 0 \f]
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 * @note
 * The Mass, Derivative and LIFT matrix for each order is include from the header files.
 * In the future, these matrix will be generate from library.
 *
 * @todo
 * 1. add boundary conditions
 */
int main(int argc, char **argv){

    /* read mesh */
    Mesh * mesh;
    Ncfile * outfile;
    int procid, nprocs;
    double dt, FinalTime = 2.4;
    char casename[16] = "Convection2d";

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
#ifdef TRI
        printf("\n    Tri Ele = %d \n", Ne);
#else
        printf("\n   Quad Ele = %d \n", Ne);
#endif
        printf("\n");
        printf("--------------------------------\n");
    }

#if defined TRI

    mesh = ReadTriMesh();
#elif defined QUAD

    mesh = ReadQuadMesh();
#endif

    /* redistribute element */
    LoadBalance(mesh);

    /* find element connections */
    FacePair(mesh);

    /* setup mesh */
#if defined TRI

    SetupTriCoeff(mesh);

#elif defined QUAD

    SetupQuadCoeff(mesh);
#endif

    /* set node connection */
    BuildMaps(mesh);

#if defined DEBUG /* check mesh */

    PrintMeshConnection(mesh);
#endif

    /* init mesh coeff */
    dt = InitMeshInfo(mesh, p_Nfields);
    dt = .5*dt/((p_N+1)*(p_N+1)); // CFL

    /* initial conditions */
    InitData(mesh);

    /* setup output file */
    outfile = SetupOutput(mesh, casename);

    /* solve */
    ConvectionRun2d(mesh, outfile ,FinalTime, dt);

    /* post process */
    Postprocess(mesh);

    /* finish */
    ConvectionFinish(mesh, outfile);

    return 0;
}


/**
 * @brief
 * Program Finalize
 *
 * @details
 * Deallocate all the mem and close the NetCDF file
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 * @warning
 * @attention
 * @note
 * @todo
 */
void ConvectionFinish(Mesh * mesh, Ncfile * outfile){

    int ret;
    /* close NetCDF output file */
    ret = ncmpi_close(outfile->ncfile);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    /* NcOutput.c */
    free(outfile->varid);

    /* Mesh2d.c */
    DestroyMatrix(mesh->GX);
    DestroyMatrix(mesh->GY);
    DestroyIntMatrix(mesh->EToV);

    /* PairFace */
    DestroyIntVector(mesh->Npar);
    DestroyIntMatrix(mesh->EToE);
    DestroyIntMatrix(mesh->EToF);
    DestroyIntMatrix(mesh->EToP);
    free(mesh->parK);
    free(mesh->parF);

    /* SetUp.c */
    DestroyVector(mesh->r);
    DestroyVector(mesh->s);
    DestroyVector(mesh->w);
    DestroyVector(mesh->wv);
    DestroyMatrix(mesh->Dr);
    DestroyMatrix(mesh->Ds);
    DestroyMatrix(mesh->LIFT);
    DestroyIntMatrix(mesh->Fmask);
    DestroyVector(mesh->rk4a);
    DestroyVector(mesh->rk4b);
    DestroyVector(mesh->rk4c);

    /* BuildMaps.c */
    DestroyIntVector(mesh->vmapM);
    DestroyIntVector(mesh->vmapP);
    DestroyIntVector(mesh->parmapOUT);
    free(mesh->f_outQ);
    free(mesh->f_inQ);

    /* InitialCondition.c */
    free(mesh->f_Q);
    free(mesh->f_rhsQ);
    free(mesh->f_resQ);
    free(mesh->f_s);

    /* InitMeshInfo.c */
    free(mesh->f_LIFT);
    free(mesh->f_Dr);
    free(mesh->f_Ds);
    free(mesh->vgeo);
    free(mesh->surfinfo);
    free(mesh->ciradius);


    MPI_Finalize();
}