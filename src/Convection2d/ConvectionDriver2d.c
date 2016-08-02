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
 *     mpirun -n 2 ./Convection2d 2 80
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 * @todo
 * 1. Add slope limiter
 * 2. Add boundary conditions
 */

#include "ConvectionDriver2d.h"

void str2int(char *str, int *N, char* errmessage);

int main(int argc, char **argv){
    int N, Ne;     /* parameters */

    double dt, FinalTime = 2.4;         /* time */
    char casename[16] = "Convection2d"; /* output filename */

    int procid, nprocs; /* process number */

    MPI_Init(&argc, &argv);                 /* initialize MPI */
    MPI_Comm_rank(MPI_COMM_WORLD, &procid); /* read process id */
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); /* read process num */

    /* get degree and element number */
    str2int(argv[1], &N, "Wrong degree input");
    str2int(argv[2], &Ne, "Wrong element number input");

    if(!procid) {
        printf("--------------------------------\n");
        printf("          Convection2d\n");
        printf("--------------------------------\n");
        printf("\n    2d Convection Test Case\n");
        printf("\n        Deg = %d \n", N);
        printf("\n    Tri Ele = %d \n", Ne);
        printf("\n   Quad Ele = %d \n", Ne);
        printf("\n");
        printf("--------------------------------\n");
    }

    StdRegions2d *shape = GenStdTriEle(N);
//    StdRegions2d *shape = GenStdQuadEle(N);
    MultiReg2d *mesh    = ReadMesh(shape, Ne);


    /* physics */
    PhysDomain2d *phys     = GenPhysDomain2d(mesh, 1);
    PhysDomain2d *flowRate = GenPhysDomain2d(mesh, 2);

    /* init phys */
    dt = InitCondition(phys, flowRate);

#if 0
    char *sname = "scalar";
    PrintPhys(phys, sname);
    char *fname = "flowRate";
    PrintPhys(flowRate, fname)
    if(!procid) {
        printf("\ntime step = %f\n", dt);
    }
#endif

    Ncfile * outfile;
    /* setup output file */
    outfile = SetupOutput(mesh, casename);

    /* solve */
    ConvectionRun2d(phys, flowRate, outfile ,FinalTime, dt);

    /* post process */
    Postprocess(phys);

    /* finish */
    int ret;
    ret = ncmpi_close(outfile->ncfile);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    free(outfile->varid);
    FreeStdRegions2d(shape);
    FreeMultiReg2d(mesh);
    FreePhysDomain2d(phys);
    FreePhysDomain2d(flowRate);

    MPI_Finalize();

    return 0;
}

void str2int(char *str, int *N, char* errmessage){
    int info = sscanf(str,"%d",N);
    if (info!=1) {
        fprintf(stderr, "%s:%s \n", errmessage, str);
        exit(-1);
    }
}