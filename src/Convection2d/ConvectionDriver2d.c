/**
 * @brief
 * Two dimensional convection problem
 *
 * @details
 * Two dimensional scalar convection problem
 * \f[\frac{\partial C}{\partial t}+\frac{\partial uC}{\partial x}+\frac{\partial vC}{\partial y}=0\f]
 *
 * Usages:
 * Use the 2 order basis with uniform mesh of 80 elements on each edge:
 *
 *     mpirun -n 2 ./Convection2d tri 2 80
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 */

#include "ConvectionDriver2d.h"

int main(int argc, char **argv){
    int N, Ne;     /* parameters */

    double dt, FinalTime = 2.4;         /* time */
    char casename[16] = "Convection2d"; /* output filename */

    int procid, nprocs; /* process number */

    MPI_Init(&argc, &argv);                 /* initialize MPI */
    MPI_Comm_rank(MPI_COMM_WORLD, &procid); /* read process id */
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); /* read process num */

    /* get degree and element number */
    str2int(argv[2], &N, "Wrong degree input");
    str2int(argv[3], &Ne, "Wrong element number input");

    if(!procid) {
        printf("--------------------------------\n");
        printf("          Convection2d\n");
        printf("--------------------------------\n");
        printf("\n    2d Convection Test Case\n");
        printf("\n        Deg = %d \n", N);
        printf("\n    Tri Ele = %d \n", Ne);
        printf("\n   Quad Ele = %d \n", Ne);
        printf("\n   Ele Type = %s \n", argv[1]);
        printf("\n");
        printf("--------------------------------\n");
    }

    StdRegions2d *shape;

    if ( !(memcmp(argv[1], "tri", 3)) ){
        shape = GenStdTriEle(N);
    }else if( !(memcmp(argv[1], "quad", 4)) ){
        shape = GenStdQuadEle(N);
    }else{
        printf("Wrong mesh type: %s.\n"
                       "The input should be either \'tri\' or \'quad\'.\n", argv[1]);
        MPI_Finalize(); exit(1);
    }

    /* gen unifrom unsructed mesh */
    MultiReg2d *mesh    = ReadMesh(shape, Ne);

    /* physics */
    PhysDomain2d *phys     = GenPhysDomain2d(mesh, 1);
    PhysDomain2d *flowRate = GenPhysDomain2d(mesh, 2);

    /* init phys */
    dt = InitCondition(phys, flowRate);

    /* setup output file */
    NcFile *outfile = SetupOutput(mesh, casename);

    /* solve */
    ConvectionRun2d(phys, flowRate, outfile ,FinalTime, dt);

    /* post process */
    Postprocess(phys);

    /* finish */
    FreeNcFile(outfile);
    FreeStdRegions2d(shape);
    FreeMultiReg2d(mesh);
    FreePhysDomain2d(phys);
    FreePhysDomain2d(flowRate);

    MPI_Finalize();

    return 0;
}