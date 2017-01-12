//
// Created by li12242 on 17/1/12.
//

#include "swe_getparameter.h"
#include "LibUtilities/UTest.h"
#include "swe_dirver2d.h"

char filename[] = "swe2d_parameter.inc";
int linenum_problem, linenum_sc, linenum_mesh;
int linenum_time, linenum_phys, linenum_LDG;

static void swe_writeInputFile();
static void swe_readInputFile();
void writeBlankine(FILE *file, int linnum);

/**
 * @brief read command line parameters
 * @param [in] argc No. of input parameter
 * @param [in] argv character pointers of input parameter
 */
void swe_getparameter(int argc, char **argv){

    int ishelp, isverbose;
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* read from command */
    UTest_Command(argc, argv, &ishelp, &isverbose);

    /* const strings */
    char helpinfo[] = HEADFINISH "DGOM:\n" HEADLINE "2d swe solver\n"
            HEADLINE "Optional features:\n"
            HEADLINE "   -help         print help information\n"
            HEADLINE "   -preprocess   generate input file\n"
            HEADLINE "Example usages:\n"
            HEADLINE "   mpirun -n 2 -host localhost ./swe2d -preprocess\n"
            HEADFINISH "\n";

    if(ishelp & (!procid)){
        printf("%s", helpinfo); exit(0);
    }

    /* generate standard input file */
    int i;
    for(i=0;i<argc;i++){
        if(!(memcmp(argv[i], "-preprocess", 11)) ){
            if(!procid){ swe_writeInputFile(); }
            exit(0);
        }
    }

    /* read input file */
    swe_readInputFile();
    return;
}


/**
 * @brief read the parameter from input file
 */
static void swe_readInputFile(){

    extern swe_solver solver;
    FILE *fp = fopen(filename, "r");

    const int len = 200;
    char buffer[len];
    int i, id, linenum=0;
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // read problem info
    for(i=0;i<linenum_problem;i++){ fgets(buffer, len, fp); linenum++; }
    // read problem indicator
    fscanf(fp, "%d\n", &id); linenum++; solver.caseid = id;

    // check element type
    if( (solver.caseid!=swe_dambreakdry) & (solver.caseid!=swe_dambreakwet)
        & (solver.caseid!=swe_parabolicbowl) & (solver.caseid!=swe_userset) ){

        if(!procid){
            fprintf(stderr, HEADFINISH "%s:\n"
                            HEADLINE " Line %d: case id %d fault. \n"
                            HEADLINE " The input type indicator should be one of \n"
                            HEADLINE "   %d - swe_dambreakdry\n"
                            HEADLINE "   %d - swe_dambreakwet\n"
                            HEADLINE "   %d - swe_parabolicbowl\n"
                            HEADLINE "   %d - swe_userset\n",
                    __FILE__, linenum, id,
                    swe_dambreakdry, swe_dambreakwet,
                    swe_parabolicbowl, swe_userset);
        }
        exit(-1);
    }

    // read cell info
    for(i=0;i<linenum_sc;i++){ fgets(buffer, len, fp); linenum++; }
    fscanf(fp, "%d\n", &id); linenum++; solver.celltype = (sc_cellType) id;

    // check element type
    if( (solver.celltype!=TRIANGLE) && (solver.celltype!=QUADRIL) ){
        fprintf(stderr, HEADFINISH "%s:\n"
                        HEADLINE " Line %d: element type %d fault. \n"
                        HEADLINE " The input type indicator should be one of \n"
                        HEADLINE "   %d - tri\n"
                        HEADLINE "   %d - quad\n",
                __FILE__, linenum, id, TRIANGLE, QUADRIL);
        exit(-1);
    }
    fscanf(fp, "%d\n", &id); linenum++; solver.N = id;

    // read mesh info
    for(i=0;i<linenum_mesh;i++){ fgets(buffer, len, fp); linenum++; }
    fscanf(fp, "%d\n", &id); linenum++; solver.Mx = id;
    fscanf(fp, "%d\n", &id); linenum++; solver.My = id;
    fscanf(fp, "%s\n", buffer); linenum++;
    solver.meshname = (char *)malloc(strlen(buffer)*sizeof(char));
    strcpy(solver.meshname, buffer);

    double val;
    // read time info
    for(i=0;i<linenum_time;i++){ fgets(buffer, len, fp); linenum++; }
    fscanf(fp, "%lf\n", &val); linenum++; solver.cfl = val;
#if DEBUG
    printf("cfl = %f\n", val);
#endif

    if( solver.cfl<0.0 | solver.cfl>1.0 ){
        fprintf(stderr, HEADFINISH "%s:\n"
                        HEADLINE " Line %d: CFL number %f fault.\n"
                        HEADLINE " The input CFL number should between [0,1] \n",
                __FILE__, linenum, val);
        exit(-1);
    }
    fscanf(fp, "%lf\n", &val); linenum++; solver.ftime = val;
    if( solver.ftime<0.0 ){
        fprintf(stderr, HEADFINISH "%s:\n"
                        HEADLINE " Line %d: final time %f fault.\n"
                        HEADLINE " The input Final time should be positive\n",
                __FILE__, linenum, val);
        exit(-1);
    }

    // read phys info
    for(i=0;i<linenum_phys;i++){ fgets(buffer, len, fp); linenum++; }
    fscanf(fp, "%lf\n", &val); linenum++; solver.gra = val;
    fscanf(fp, "%lf\n", &val); linenum++; solver.hcrit = val;
    fscanf(fp, "%lf\n", &val); linenum++; solver.roughness = val;

    // real LDG parameters
    for(i=0;i<3;i++){ fgets(buffer, len, fp); linenum++; }
    fscanf(fp, "%lf\n", &val); linenum++; solver.LDG_parameter[0] = val;
    fscanf(fp, "%lf\n", &val); linenum++; solver.LDG_parameter[1] = val;
    fscanf(fp, "%lf\n", &val); linenum++; solver.LDG_parameter[2] = val;
}

/**
 * @brief generate standard input file
 */
static void swe_writeInputFile(){

#undef HEADLINE
#undef HEADFINISH

#define HEADFINISH "[============]"
#define HEADLINE   "[------------]"

    char probleminfo[] = HEADFINISH "DGOM: 2d swe solver\n"
            HEADLINE "case indicator (1 parameter)\n"
            HEADLINE "    1. case indicator  |-- 1. dam break wet\n"
            HEADLINE "                       |-- 2. dam break dry\n";
            HEADLINE "                       |-- 3. parabolic bowl\n";
            HEADLINE "                       |-- 4. user specific\n";

    char scinfo[] = HEADFINISH "standard element info (2 parameters)\n"
            HEADLINE "    1. element type  |--- 0. triangle\n"
            HEADLINE "                     |--- 1. quadrilateral\n"
            HEADLINE "    2. order of polynomial\n";

    char meshinfo[] = HEADFINISH "mesh info (2 parameters)\n"
            HEADLINE "    1. num of elements in x direction\n"
            HEADLINE "    2. num of elements in y direction\n";
            HEADLINE "    3. user specific mesh name\n";

    char timeinfo[] = HEADFINISH "time info (2 parameters)\n"
            HEADLINE "    1. CFL number\n"
            HEADLINE "    2. final time\n";

    char physinfo[] = HEADFINISH "physical parameter (3 parameters)\n"
            HEADLINE "    1. gravity acceleration\n"
            HEADLINE "    2. hcrit -- minimum water depth of wet cell\n"
            HEADLINE "    3. manning coefficient for friction\n";
            HEADLINE "    4. file name of topography elevation\n";

    char LDGinfo[] = HEADFINISH "LDG parameter for viscosity term (3 parameters)\n"
            HEADLINE "    1. C11\n"
            HEADLINE "    2. C12\n"
            HEADLINE "    3. C22\n";

    FILE *wfile = fopen(filename, "w");

    fprintf(wfile, "%s", probleminfo); linenum_problem = 6;
    writeBlankine(wfile, 1);

    fprintf(wfile, "%s", scinfo); linenum_sc = 4;
    writeBlankine(wfile, 2);

    fprintf(wfile, "%s", meshinfo); linenum_mesh = 4;
    writeBlankine(wfile, 3);

    fprintf(wfile, "%s", timeinfo); linenum_time = 3;
    writeBlankine(wfile, 2);

    fprintf(wfile, "%s", physinfo); linenum_phys = 5;
    writeBlankine(wfile, 4);

    fprintf(wfile, "%s", LDGinfo); linenum_LDG = 4;
    writeBlankine(wfile, 3);

    fclose(wfile);
}

/**
 * @brief write blank line to file
 */
static void writeBlankine(FILE *file, int linnum){
    register int i;
    for(i=0;i<linnum;i++)
        fprintf(file, "\n");
}
