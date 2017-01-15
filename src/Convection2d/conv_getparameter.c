//
// Created by li12242 on 16/12/27.
//

#include "conv_getparameter.h"
#include "LibUtilities/UTest.h"
#include "conv_driver2d.h"

#define DEBUG 0

char FILE_NAME[] = "conv2d_paramter.inc";

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

    /* const strings */
    char helpinfo[] = HEADEND "DGOM:\n" HEADLINE "2d convection problem\n"
            HEADLINE "Optional features:\n"
            HEADLINE "   -help         print help information\n"
            HEADLINE "   -preprocess   generate input file\n"
            HEADLINE "Example usages:\n"
            HEADLINE "   mpirun -n 2 -host localhost ./convection2d -preprocess\n"
            HEADEND "\n";

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
    FILE *fp = fopen(FILE_NAME, "r");

    const int len = 200;
    char buffer[len];
    int i, integer, linenum=0;
    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    // read problem info
    for(i=0;i<4;i++){ fgets(buffer, len, fp); linenum++; }
    // read problem indicator
    fscanf(fp, "%d\n", &integer); linenum++;
    solver.caseid = integer;

    if(!procid){
        printf(HEADEND "--------------------------------\n");
        printf(HEADLINE "          Convection2d\n");
        printf(HEADLINE "--------------------------------\n");
        switch (integer){
            case conv_rotational_convection:
                printf(HEADLINE " case: rotational convection\n"); break;
            case conv_advection_diffusion:
                printf(HEADLINE " case: advection diffusion\n"); break;
            default:
                fprintf(stderr, HEADEND "%s:\n"
                                HEADLINE " Line %d: case id %d fault. \n"
                                HEADLINE " The input type indicator should be one of \n"
                                HEADLINE "   %d - rotational convection\n"
                                HEADLINE "   %d - advection-diffusion\n",
                        __FILE__, linenum, integer,
                        conv_rotational_convection,
                        conv_advection_diffusion);
                MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

    // read element info
    for(i=0;i<4;i++){ fgets(buffer, len, fp); linenum++; }

    fscanf(fp, "%d\n", &integer); linenum++;
    solver.celltype = (sc_cellType) integer;
    fscanf(fp, "%d\n", &integer); linenum++;
    solver.N = integer;

    if(!procid){
        switch (solver.celltype){
            case TRIANGLE:
                printf(HEADLINE " cell type: triangle\n"); break;
            case QUADRIL:
                printf(HEADLINE " cell type: quadrilateral\n"); break;
            default:
                fprintf(stderr, HEADEND "%s:\n"
                                HEADLINE " Line %d: element type %d fault. \n"
                                HEADLINE " The input type indicator should be one of \n"
                                HEADLINE "   %d - tri\n"
                                HEADLINE "   %d - quad\n",
                        __FILE__, linenum, integer, TRIANGLE, QUADRIL);
                MPI_Abort(MPI_COMM_WORLD, -1);
        }
        printf(HEADLINE " polynomial degree: %d\n", solver.N);
    }

    // read mesh info
    for(i=0;i<3;i++){ fgets(buffer, len, fp); linenum++; }
    fscanf(fp, "%d\n", &integer); linenum++;
    fscanf(fp, "%d\n", &integer); linenum++;
    solver.Ne = integer;

    if(!procid){
        printf(HEADLINE " number of cell on x: %d\n", solver.Ne);
        printf(HEADLINE " number of cell on y: %d\n", solver.Ne);
    }

    float val;
    // read time info
    for(i=0;i<3;i++){ fgets(buffer, len, fp); linenum++; }
    fscanf(fp, "%f\n", &val); linenum++;
    solver.cfl = val;

    if(!procid){
        if( solver.cfl>0.0 | solver.cfl<1.0 ){
            printf(HEADLINE " cfl = %f\n", solver.cfl);
        }else{
            fprintf(stderr, HEADEND "%s:\n"
                            HEADLINE " Line %d: CFL number %f fault.\n"
                            HEADLINE " The input CFL number should between [0,1] \n",
                    __FILE__, linenum, val);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

    fscanf(fp, "%f\n", &val); linenum++;
    solver.finaltime = val;

    if(!procid){
        // check final time
        if( val>0.0 ){
            printf(HEADLINE " final time = %f\n", solver.finaltime);
        }else{
            fprintf(stderr, HEADEND "%s:\n"
                            HEADLINE " Line %d: final time %f fault.\n"
                            HEADLINE " The input Final time should be positive\n",
                    __FILE__, linenum, val);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

    // read phys info
    for(i=0;i<4;i++){ fgets(buffer, len, fp); linenum++; }
    // read phys parameter
    fscanf(fp, "%f\n", &val); linenum++; solver.u = val;
    fscanf(fp, "%f\n", &val); linenum++; solver.v = val;
    fscanf(fp, "%f\n", &val); linenum++; solver.viscosity = val;

    if(!procid){
        printf(HEADLINE " u = %f\n", solver.u);
        printf(HEADLINE " v = %f\n", solver.v);
        printf(HEADLINE " vis = %f\n", solver.viscosity);
    }

    // real LDG parameters
    for(i=0;i<4;i++){ fgets(buffer, len, fp); linenum++; }
    fscanf(fp, "%f\n", &val); linenum++; solver.LDG_parameter[0] = val;
    fscanf(fp, "%f\n", &val); linenum++; solver.LDG_parameter[1] = val;
    fscanf(fp, "%f\n", &val); linenum++; solver.LDG_parameter[2] = val;

    if(!procid){
        printf(HEADLINE " c11 = %f\n", solver.LDG_parameter[0]);
        printf(HEADLINE " c12 = %f\n", solver.LDG_parameter[1]);
        printf(HEADLINE " c22 = %f\n", solver.LDG_parameter[2]);
    }
    fclose(fp);
}

/**
 * @brief generate the input file
 */
static void conv_genInputFile(){

#undef HEADLINE
#undef HEADEND

#define HEADEND   "[============]"
#define HEADLINE  "[------------]"

    char probleminfo[] = HEADEND "DGOM: 2d convection problem\n"
            HEADLINE "case indicator (1 parameter)\n"
            HEADLINE "    1. case indicator  |-- 1. rotational convection\n"
            HEADLINE "                       |-- 2. advection-diffusion\n";

    char scinfo[] = HEADEND "standard element info (2 parameters)\n"
            HEADLINE "    1. element type  |--- 0. triangle\n"
            HEADLINE "                     |--- 1. quadrilateral\n"
            HEADLINE "    2. order of polynomial\n";

    char meshinfo[] = HEADEND "mesh info (2 parameters)\n"
            HEADLINE "    1. num of elements in x direction\n"
            HEADLINE "    2. num of elements in y direction\n";

    char timeinfo[] = HEADEND "time info (2 parameters)\n"
            HEADLINE "    1. CFL number\n"
            HEADLINE "    2. final time\n";

    char physinfo[] = HEADEND "physical parameter for advection-diffusion case (3 parameters)\n"
            HEADLINE "    1. flow rate u for x direction\n"
            HEADLINE "    2. flow rate v for x direction\n"
            HEADLINE "    3. viscosity parameter miu\n";

    char LDGinfo[] = HEADEND "LDG parameter for viscosity term (3 parameters)\n"
            HEADLINE "    1. C11 (0.0 recommend)\n"
            HEADLINE "    2. C12 (0.5 recommend)\n"
            HEADLINE "    3. C22\n";

    FILE *wfile = fopen(FILE_NAME, "w");

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

    fprintf(wfile, "%s", physinfo);
    fprintf(wfile, "\n");
    fprintf(wfile, "\n");
    fprintf(wfile, "\n");

    fprintf(wfile, "%s", LDGinfo);
    fprintf(wfile, "\n");
    fprintf(wfile, "\n");
    fprintf(wfile, "\n");

    fclose(wfile);
}

