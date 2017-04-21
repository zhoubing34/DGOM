//
// Created by li12242 on 17/1/12.
//

#include "../SWELib/swe_lib.h"

#define DEBUG 0

static SWE_Run_Type swe_read_command(int argc, char **argv);
static void swe_create_inputfile(char *filename);

void swe_input(int argc, char **argv){
    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    SWE_Run_Type run_type = swe_read_command(argc, argv);
    char helpinfo[] = HEAD_LINE "DGOM:\n" HEAD_LINE "2d swe problem\n"
            HEAD_LINE "Optional features:\n"
            HEAD_LINE "   -help          print help information;\n"
            HEAD_LINE "   -create_input  create input file;\n"
            HEAD_LINE "   -run           run with input file;\n"
            HEAD_LINE "Example usages:\n"
            HEAD_LINE "   mpirun -n 2 -host localhost ./swe2d -create_input swe_parameter.inc\n"
            HEAD_LINE "\n"
            HEAD_LINE "   mpirun -n 2 -host localhost ./swe2d -run swe_parameter.inc\n";

    extern SWE_Solver solver;
    switch (run_type){
        case SWE_HELP:
            if(!procid) { printf("%s", helpinfo); } exit(0);
        case SWE_CREATE_INPUT:
            if(!procid) { swe_create_inputfile(argv[2]); } exit(0);
        case SWE_RUN:
            strcpy(solver.filename, argv[2]); // set the input file name
            printf(HEAD_LINE " input file: %s\n", solver.filename);
            break;
        case SWE_UNKNOWN:
            fprintf(stderr, "%s (%d), Unknown input parameter\n", __FUNCTION__, __LINE__);
    }
    return;
}

static SWE_Run_Type swe_read_command(int argc, char **argv){
    SWE_Run_Type run_type = SWE_UNKNOWN; // default

    char help_str[] = "-help";
    char pre_str[]  = "-create_input";
    char run_str[]  = "-run";

    register int i;
    for(i=0;i<argc;i++){
        if(!(memcmp(argv[i], help_str, strlen(help_str))) ){
            run_type = SWE_HELP;
        }else if(!(memcmp(argv[i], pre_str, strlen(pre_str) ))){
            run_type = SWE_CREATE_INPUT;
        }else if(!(memcmp(argv[i], run_str, strlen(run_str) ))){
            run_type = SWE_RUN;
        }
    }
    return run_type;
}

static arg_section** swe_create_section(){
    int section_num = SECTION_NUM;
    arg_section** section_p = (arg_section**) calloc((size_t)section_num, sizeof(arg_section*));

    int ind = 0;
    /// 0. section: case info
    char caseinfo[] = HEAD_LINE "DGOM: 2d swe solver\n"
            HEAD_LINE "case info (1 parameter):\n"
            HEAD_LINE "    1. case name;\n";
    int var_num = 1;
    section_p[ind++] = section_create(caseinfo, var_num);
    /// 1. section: std cell info
    char scinfo[] = HEAD_LINE "standard element info (2 parameters):\n"
            HEAD_LINE "    1. element type  |--- 2. triangle\n"
            HEAD_LINE "                     |--- 3. quadrilateral\n"
            HEAD_LINE "    2. order of polynomial\n";
    var_num = 2;
    section_p[ind++] = section_create(scinfo, var_num);
    /// 2. section: mesh info
    char meshinfo[] = HEAD_LINE "mesh info (1 parameters):\n"
            HEAD_LINE "    1. open boundary file;\n";
    var_num = 1;
    section_p[ind++] = section_create(meshinfo, var_num);
    /// 3. section: time info
    char timeinfo[] = HEAD_LINE "time info (2 parameters):\n"
            HEAD_LINE "    1. CFL number;\n"
            HEAD_LINE "    2. dt;\n"
            HEAD_LINE "    3. final time;\n";
    var_num = 3;
    section_p[ind++] = section_create(timeinfo, var_num);
    /// 4. section: physical info
    char physinfo[] = HEAD_LINE "physical parameter (3 parameters)\n"
            HEAD_LINE "    1. gravity acceleration;\n"
            HEAD_LINE "    2. hcrit -- minimum water depth of wet cell;\n"
            HEAD_LINE "    3. manning coefficient file;\n"
            HEAD_LINE "    4. initial condition file (h, hu, hv, b);\n";
    var_num = 4;
    section_p[ind++] = section_create(physinfo, var_num);
    /// 5. section: result info
    char outinfo[] = HEAD_LINE "information for output (2 parameters)\n"
            HEAD_LINE "    1. filename;\n"
            HEAD_LINE "    2. output time interval;\n";
    var_num = 2;
    section_p[ind++] = section_create(outinfo, var_num);

    return section_p;
}

void swe_free_section(arg_section **section_p){
    int n;
    for(n=0;n<SECTION_NUM;n++) {section_free(section_p[n]);}
    free(section_p);
    return;
}

/**
 * @brief create input file for swe2d solver.
 * @param filename name of input file;
 */
static void swe_create_inputfile(char *filename){
    arg_section **section_p = swe_create_section();

    int n;
    for(n=0;n<SECTION_NUM;n++) {section_print(section_p[n]);}

    FILE *fp = fopen(filename, "w");
    for(n=0;n<SECTION_NUM;n++){
        section_write_file(fp, section_p[n]);
    }
    fclose(fp);
    swe_free_section(section_p);
}

arg_section** swe_read_section(){
    arg_section **arg = swe_create_section();
    extern SWE_Solver solver;
    FILE *fp = fopen(solver.filename, "r");
    int n;
    for(n=0;n<SECTION_NUM;n++){ section_read_file(fp, arg[n]); }
    fclose(fp);
    return arg;
}

//SWE_Solver* swe_create_solver(){
//
//    arg_section **arg = swe_read_section();
//    SWE_Solver *solver = (SWE_Solver*) calloc(1, sizeof(SWE_Solver));
//
//    int procid;
//    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
//
//    /// 0. section: case info
//    arg_section *sec_p = arg[0];
//    sscanf(sec_p->arg_vec_p[0], "%d\n", &(solver->caseid));
//
//    if(!procid){
//        if( (solver->caseid==swe_dambreakdry) || (solver->caseid==swe_dambreakwet)
//            || (solver->caseid==swe_parabolicbowl) || (solver->caseid==swe_userset) ){
//            printf(HEAD_TITLE " DGOM: 2d swe solver\n"
//                           HEAD_LINE " case = %d\n", solver->caseid);
//        }else{
//            fprintf(stderr, "%s (%d): Unknown case id %d\n",
//                    __FILE__, __LINE__, solver->caseid);
//            MPI_Abort(MPI_COMM_WORLD, -1);
//        }
//    }
//    /// 1. section: cell info
//    sec_p = arg[1];
//    sscanf(sec_p->arg_vec_p[0], "%d\n", &(solver->celltype));
//    sscanf(sec_p->arg_vec_p[1], "%d\n", &(solver->N));
//
//    if(!procid){
//        if( solver->celltype==TRIANGLE ){
//            printf(HEAD_LINE " cell type: triangle\n");
//        }else if(solver->celltype==QUADRIL){
//            printf(HEAD_LINE " cell type: quadrilateral\n");
//        }else{
//            fprintf(stderr, "%s (%d): Unknown cell type %d\n",
//                    __FILE__, __LINE__, solver->celltype);
//            MPI_Abort(MPI_COMM_WORLD, -1);
//        }
//    }
//    /// 2. section: mesh info
//    sec_p = arg[2];
//    int Mx, My;
//    double xmin, xmax, ymin, ymax;
//    sscanf(sec_p->arg_vec_p[0], "%d\n", &(Mx));
//    sscanf(sec_p->arg_vec_p[1], "%d\n", &(My));
//    solver->casename = calloc(strlen(sec_p->arg_vec_p[2])+1, sizeof(char));
//    strcpy(solver->casename, sec_p->arg_vec_p[2]);
//    char *casename = solver->casename;
//    switch (solver->caseid){
//        case swe_dambreakwet:
//            xmin=0, xmax=1000, ymin=-100, ymax=100;
//            if(!procid){ printf(HEAD_LINE " cell num: Mx=%d, My=%d\n", Mx, My); }
//            solver->phys = swe_uniform_mesh(solver, Mx, My, xmin, xmax, ymin, ymax);
//            break;
//        case swe_dambreakdry:
//            xmin=0, xmax=1000, ymin=-100, ymax=100;
//            if(!procid){ printf(HEAD_LINE " cell num: Mx=%d, My=%d\n", Mx, My); }
//            solver->phys = swe_uniform_mesh(solver, Mx, My, xmin, xmax, ymin, ymax);
//            break;
//        case swe_parabolicbowl:
//            xmin = -4000, xmax = 4000, ymin = -4000, ymax = 4000;
//            if(!procid){  printf(HEAD_LINE " cell num: Mx=%d, My=%d\n", Mx, My); }
//            solver->phys = swe_uniform_mesh(solver, Mx, My, xmin, xmax, ymin, ymax);
//            break;
//        case swe_userset:
//            if(!procid){  printf(HEAD_LINE " case name: %s\n", casename); }
//            solver->phys = swe_file_mesh(solver, casename);
//            break;
//        default:
//            MPI_Abort(MPI_COMM_WORLD, -1);
//    }
//    /// 3. section: time info
//    sec_p = arg[3];
//    sscanf(sec_p->arg_vec_p[0], "%lf\n", &(solver->cfl));
//    sscanf(sec_p->arg_vec_p[1], "%lf\n", &(solver->ftime));
//    if(!procid) printf(HEAD_LINE " cfl = %f\n", solver->cfl);
//    if(!procid) printf(HEAD_LINE " final time = %f\n", solver->ftime);
//    /// 4. section: phys info
//    sec_p = arg[4];
//    sscanf(sec_p->arg_vec_p[0], "%lf\n", &(solver->gra));
//    sscanf(sec_p->arg_vec_p[1], "%lf\n", &(solver->hcrit));
//    sscanf(sec_p->arg_vec_p[2], "%lf\n", &(solver->m));
//
//    if(!procid) printf(HEAD_LINE " gra = %f\n", solver->gra);
//    if(!procid) printf(HEAD_LINE " hcrit = %f\n", solver->hcrit);
//    if(!procid) printf(HEAD_LINE " manning factor = %f\n", solver->m);
//
//    /// 5. section: LDG info
//    sec_p = arg[5];
//    sscanf(sec_p->arg_vec_p[0], "%lf\n", &(solver->c12));
//    if(!procid) printf(HEAD_LINE " c12 = %f\n", solver->c12);
//
//    swe_free_section(arg);
//    return solver;
//}
