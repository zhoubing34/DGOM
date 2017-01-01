//
// Created by li12242 on 16/12/27.
//

#include "conv_getparameter.h"
#include "LibUtilities/UTest.h"
#include "conv_driver2d.h"

void conv_getparameter(int argc, char **argv){

    int ishelp, isverbose;

    UTest_Command(argc, argv, &ishelp, &isverbose);

    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if(ishelp & (!procid)){
        printf(HEADFINISH "DGOM:\n" HEADLINE "2d convection problem\n"
                       HEADLINE "Command parameters: \n"
                       HEADLINE "   N - order \n"
                       HEADLINE "   Ne - cell number\n"
                       HEADLINE "   Cell Type - tri=0/quad=1\n"
                       HEADLINE "   CFL - CFL number\n"
                       HEADLINE "   time - final time\n"
                       HEADLINE "Example usages:\n"
                       HEADLINE "   mpirun -n 2 -host localhost ./convection2d  2  80  0  0.3  2.4\n"
                       HEADLINE "\n"
                       HEADLINE "Optional features:\n"
                       HEADLINE "   -help     print help information\n"
                       HEADFINISH "   -verbose  print variables to log files\n\n");
        exit(0);
    }

    /* save to global variable */
    extern conv_solver2d solver;
    solver.isverbose = isverbose;

    int N, Ne, type, par=1;
    double CFL,ftime;
    str2int(argv[par++], &N, "order");
    str2int(argv[par++], &Ne, "number of cells");
    str2int(argv[par++], &type, "cell type");
    str2double(argv[par++], &CFL, "CFL number");
    str2double(argv[par++], &ftime, "final time");

    /* check parameters */
    if( (type!=0) && (type!=1) ){
        fprintf(stderr, HEADFINISH "%s:\n"
                       HEADLINE " Element type %d fault. \n"
                       HEADLINE " The input type indicator should be one of \n"
                       HEADLINE "   0 - tri\n"
                       HEADLINE "   1 - quad\n", __FILE__, type);
        exit(-1);
    }

    if( CFL<0.0 | CFL>1.0 ){
        fprintf(stderr, HEADFINISH "%s:\n"
                       HEADLINE " CFL number %f fault.\n "
                       HEADLINE " The input CFL number should between [0,1] \n", __FILE__, CFL);
        exit(-1);
    }

    if( ftime<0.0 ){
        fprintf(stderr, HEADFINISH "%s:\n"
                       HEADLINE " Final time %f fault.\n "
                       HEADLINE " The input Final time should be positive\n", __FILE__, ftime);
        exit(-1);
    }


    /* assignment */
    solver.N = N;
    solver.Ne = Ne;
    solver.type = (sc_cellType) type;
    solver.cfl = CFL;
    solver.finalTime = ftime;

    return;
}

