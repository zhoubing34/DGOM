//
// Created by li12242 on 16/12/27.
//

#include "conv_getparameter.h"
#include "LibUtilities/UTest.h"
#include "conv_driver2d.h"

#define DEBUG 0

/* const strings */
char helpinfo[] = HEADFINISH "DGOM:\n" HEADLINE "2d convection problem\n"
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
        HEADFINISH "   -verbose  print variables to log files\n\n";

char filename[] = "conv2d_paramter.inc";

/* generate the input file */
static void conv_genInputFile();
/* read input file */
static void conv_readInputFile();

/**
 * @brief read parameters from input command or file
 * @param [in] argc
 * @param [in] argv
 */
void conv_getparameter(int argc, char **argv){

    int ishelp, isverbose;

    /* read from command */
    UTest_Command(argc, argv, &ishelp, &isverbose);

    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if(ishelp & (!procid)){
        printf("%s", helpinfo);
        exit(0);
    }

    /* save to global variable */
    extern conv_solver2d solver;
    solver.isverbose = isverbose;

    /* generate standard input file */
    int i;
    for(i=0;i<argc;i++){
        if(!(memcmp(argv[i], "-preprocess", 11)) ){
            if(!procid){ conv_genInputFile(); }
            exit(0);
        }
    }

    /* read input file */
    conv_readInputFile();
    return;
}

/**
 * @brief read parameters from input file
 */
static void conv_readInputFile(){

    /* save to global variable */
    extern conv_solver2d solver;
    FILE *fp = fopen(filename, "r");

    const int len = 200;
    char buffer[len];
    int i, integer, linenum=0;

    // read problem info
    for(i=0;i<4;i++){ fgets(buffer, len, fp); linenum++; }
    // read problem indicator
    fscanf(fp, "%d\n", &integer); linenum++;
    solver.caseid = integer;

    // check element type
    if( (integer!=conv_rotational_convection) && (integer!=conv_advection_diffusion) ){
        fprintf(stderr, HEADFINISH "%s:\n"
                HEADLINE " Line %d: case id %d fault. \n"
                HEADLINE " The input type indicator should be one of \n"
                HEADLINE "   %d - rotational convection\n"
                HEADLINE "   %d - advection-diffusion\n",
                __FILE__, linenum, integer,
                conv_rotational_convection,
                conv_advection_diffusion);
        exit(-1);
    }

    // read element info
    for(i=0;i<4;i++){ fgets(buffer, len, fp); linenum++; }

    fscanf(fp, "%d\n", &integer); linenum++;

    // check element type
    if( (integer!=TRIANGLE) && (integer!=QUADRIL) ){
        fprintf(stderr, HEADFINISH "%s:\n"
                HEADLINE " Line %d: element type %d fault. \n"
                HEADLINE " The input type indicator should be one of \n"
                HEADLINE "   %d - tri\n"
                HEADLINE "   %d - quad\n",
                __FILE__, linenum, integer, TRIANGLE, QUADRIL);
        exit(-1);
    }

    solver.celltype = (sc_cellType) integer;
    fscanf(fp, "%d\n", &integer); linenum++;
    solver.N = integer;

    // read mesh info
    for(i=0;i<3;i++){ fgets(buffer, len, fp); linenum++; }
    fscanf(fp, "%d\n", &integer); linenum++;
    fscanf(fp, "%d\n", &integer); linenum++;
    solver.Ne = integer;

    float val;
    // read time info
    for(i=0;i<3;i++){ fgets(buffer, len, fp); linenum++; }
    fscanf(fp, "%f\n", &val); linenum++;

#if DEBUG
    printf("cfl = %f\n", val);
#endif

    // check cfl number
    if( val<0.0 | val>1.0 ){
        fprintf(stderr, HEADFINISH "%s:\n"
                HEADLINE " Line %d: CFL number %f fault.\n"
                HEADLINE " The input CFL number should between [0,1] \n",
                __FILE__, linenum, val);
        exit(-1);
    }


    solver.cfl = val;
    fscanf(fp, "%f\n", &val); linenum++;

    // check final time
    if( val<0.0 ){
        fprintf(stderr, HEADFINISH "%s:\n"
                HEADLINE " Line %d: final time %f fault.\n"
                HEADLINE " The input Final time should be positive\n",
                __FILE__, linenum, val);
        exit(-1);
    }

    solver.finaltime = val;
    fclose(fp);
}

/**
 * @brief generate the input file
 */
static void conv_genInputFile(){

#undef HEADLINE
#undef HEADFINISH

#define HEADFINISH "[============]"
#define HEADLINE   "[------------]"

    char probleminfo[] = HEADFINISH "DGOM: 2d convection problem\n"
            HEADLINE "case indicator (1 parameter)\n"
            HEADLINE "    1. case indicator  |-- 1. rotational convection\n"
            HEADLINE "                       |-- 2. advection-diffusion\n";

    char scinfo[] = HEADFINISH "standard element info (2 parameters)\n"
            HEADLINE "    1. element type  |--- 0. triangle\n"
            HEADLINE "                     |--- 1. quadrilateral\n"
            HEADLINE "    2. order of polynomial\n";

    char meshinfo[] = HEADFINISH "mesh info (2 parameters)\n"
            HEADLINE "    1. num of elements in x direction\n"
            HEADLINE "    2. num of elements in y direction\n";

    char timeinfo[] = HEADFINISH "time info (2 parameters)\n"
            HEADLINE "    1. CFL number\n"
            HEADLINE "    2. final time\n";

    FILE *wfile = fopen(filename, "w");

    fprintf(wfile, "%s", probleminfo);
    fprintf(wfile, "\n");

    fprintf(wfile, "%s", scinfo);
    fprintf(wfile, "\n");
    fprintf(wfile, "\n");

    fprintf(wfile, "%s", meshinfo);
    fprintf(wfile, "\n");
    fprintf(wfile, "\n");

    fprintf(wfile, "%s", timeinfo);
    fprintf(wfile, "\n");
    fprintf(wfile, "\n");

    fclose(wfile);
}

