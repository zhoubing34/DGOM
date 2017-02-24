//
// Created by li12242 on 16/12/27.
//

#include "conv_driver2d.h"
#include "Utility/arg_section.h"

#define HEADEND   "[============]"
#define HEADLINE  "[------------]"
#define SEC_NUM 6

typedef enum {
    conv_command_print_help = 0,
    conv_command_create_input = 1,
    conv_command_run = 2,
}conv_command;

char FILE_NAME[] = "conv2d_paramter.inc";

/* read input file */
static void conv_read_input();
/* create input file */
static void conv_create_input();
/* read command line */
static conv_command conv_read_command(int argc, char **argv);
/* free argument section memory */
static void conv_free_section(arg_section **section_p);
/* create argument section */
static arg_section** conv_create_section();

/**
 * @brief read parameters from input command or file
 * @param [in] argc number of conmand
 * @param [in] argv pointers to the command strings
 */
void conv_getparameter(int argc, char **argv){

    /* read from command */
    conv_command run_type = conv_read_command(argc, argv);

    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    /* const strings */
    char helpinfo[] = HEADEND "DGOM:\n" HEADLINE "2d convection problem\n"
            HEADLINE "Optional features:\n"
            HEADLINE "   -help         print help information\n"
            HEADLINE "   -preprocess   generate input file\n"
            HEADLINE "Example usages:\n"
            HEADLINE "   mpirun -n 2 -host localhost ./convection2d -preprocess\n";

    if(run_type == conv_command_print_help){
        if(!procid) { printf("%s", helpinfo); } exit(0);
    }else if(run_type == conv_command_create_input){
        if(!procid) { conv_create_input(); } exit(0);
    }else if(run_type == conv_command_run){
        conv_read_input();
    }
    return;
}
/**
 * @brief create input file
 */
static void conv_create_input(){
    arg_section **sec_p = conv_create_section();
    int n;
    for(n=0;n<SEC_NUM;n++)
        section_print(sec_p[n]);

    FILE *fp = fopen(FILE_NAME, "w");
    for(n=0;n<SEC_NUM;n++){
        section_write_file(sec_p[n], fp);
    }
    fclose(fp);
    conv_free_section(sec_p);
    return;
}

/**
 * @brief read command line and get the run type.
 * @param argc number of command argument.
 * @param argv pointers to command argument.
 * @return
 * run type.
 */
static conv_command conv_read_command(int argc, char **argv){
    conv_command command_type = conv_command_run;

    char help_str[] = "-help";
    char pre_str[] = "-preprocess";

    register int i;
    for(i=0;i<argc;i++){
        if(!(memcmp(argv[i], help_str, strlen(help_str))) ){
            command_type = conv_command_print_help;
        }else if(!(memcmp(argv[i], pre_str, strlen(pre_str) ))){
            command_type = conv_command_create_input;
        }
    }
    return command_type;
}
/**
 * @brief free argument section memory.
 * @param section_p argument section strucutre.
 */
static void conv_free_section(arg_section **section_p){
    int n;
    for(n=0;n<SEC_NUM;n++)
        section_free(section_p[n]);
    free(section_p);
    return;
}
/**
 * @brief create argument section structure for convection problem.
 * @return argument section structure.
 */
static arg_section** conv_create_section(){
    arg_section** section_p = (arg_section**) calloc(SEC_NUM, sizeof(arg_section*));

    int ind = 0;
    char probleminfo[] = HEADEND "DGOM: 2d convection problem\n"
            HEADLINE "case indicator (1 parameter)\n"
            HEADLINE "    1. case indicator  |-- 1. rotational convection\n"
            HEADLINE "                       |-- 2. advection-diffusion\n";
    int var_num = 1;
    section_p[ind++] = section_create(probleminfo, var_num);

    char scinfo[] = HEADEND "standard element info (2 parameters)\n"
            HEADLINE "    1. element type  |--- 0. triangle\n"
            HEADLINE "                     |--- 1. quadrilateral\n"
            HEADLINE "    2. order of polynomial\n";
    var_num = 2;
    section_p[ind++] = section_create(scinfo, var_num);

    char meshinfo[] = HEADEND "mesh info (2 parameters)\n"
            HEADLINE "    1. num of elements in x direction\n"
            HEADLINE "    2. num of elements in y direction\n";
    var_num = 2;
    section_p[ind++] = section_create(meshinfo, var_num);

    char timeinfo[] = HEADEND "time info (2 parameters)\n"
            HEADLINE "    1. CFL number\n"
            HEADLINE "    2. dt\n"
            HEADLINE "    3. final time\n";
    var_num = 3;
    section_p[ind++] = section_create(timeinfo, var_num);

    char physinfo[] = HEADEND "physical parameter for advection-diffusion case (3 parameters)\n"
            HEADLINE "    1. flow rate u for x direction\n"
            HEADLINE "    2. flow rate v for x direction\n"
            HEADLINE "    3. viscosity parameter miu\n";
    var_num = 3;
    section_p[ind++] = section_create(physinfo, var_num);

    char LDGinfo[] = HEADEND "LDG parameter for viscosity term (3 parameters)\n"
            HEADLINE "    1. C11 (0.0 recommend)\n"
            HEADLINE "    2. C12 (0.5 recommend)\n"
            HEADLINE "    3. C22\n";
    var_num = 3;
    section_p[ind++] = section_create(LDGinfo, var_num);

    return section_p;
}
/**
 * @brief read argument section structure from input file.
 */
static void conv_read_input(){
    arg_section **sec_p = conv_create_section();
    // read section from file
    FILE *fp = fopen(FILE_NAME, "r");
    int n;
    for(n=0;n<SEC_NUM;n++){
        section_read_file(sec_p[n], fp);
    }
    fclose(fp);

    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    // global variable
    extern conv_solver2d solver;
    /// 0. section: case info
    arg_section *sec = sec_p[0];
    sscanf(sec->arg_vec_p[0], "%d\n", &(solver.caseid));

    if(!procid){
        printf(HEADEND "--------------------------------\n");
        printf(HEADLINE "          Convection2d\n");
        printf(HEADLINE "--------------------------------\n");
        switch (solver.caseid){
            case conv_rotational_convection:
                printf(HEADLINE " case: rotational convection\n"); break;
            case conv_advection_diffusion:
                printf(HEADLINE " case: advection diffusion\n"); break;
            default:
                fprintf(stderr, HEADEND "%s:\n"
                                HEADLINE " The input type indicator should be one of \n"
                                HEADLINE "   %d - rotational convection\n"
                                HEADLINE "   %d - advection-diffusion\n",
                        __FILE__, conv_rotational_convection, conv_advection_diffusion);
                MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    /// 1. cell info
    sec = sec_p[1];
    sscanf(sec->arg_vec_p[0], "%d\n", &(solver.celltype));
    sscanf(sec->arg_vec_p[1], "%d\n", &(solver.N));
    if(!procid){
        switch (solver.celltype){
            case TRIANGLE:
                printf(HEADLINE " cell type: triangle\n"); break;
            case QUADRIL:
                printf(HEADLINE " cell type: quadrilateral\n"); break;
            default:
                fprintf(stderr, HEADEND "%s:\n"
                                HEADLINE " The input type indicator should be one of \n"
                                HEADLINE "   %d - tri\n"
                                HEADLINE "   %d - quad\n",
                        __FILE__, TRIANGLE, QUADRIL);
                MPI_Abort(MPI_COMM_WORLD, -1);
        }
        printf(HEADLINE " polynomial degree: %d\n", solver.N);
    }
    /// 2. mesh info
    sec = sec_p[2];
    sscanf(sec->arg_vec_p[0], "%d\n", &(solver.Ne));
    sscanf(sec->arg_vec_p[1], "%d\n", &(solver.Ne));
    if(!procid){
        printf(HEADLINE " Ne on x: %d\n", solver.Ne);
        printf(HEADLINE " Ne on y: %d\n", solver.Ne);
    }
    /// 3. time info
    sec = sec_p[3];
    sscanf(sec->arg_vec_p[0], "%lf\n", &(solver.cfl));
    sscanf(sec->arg_vec_p[1], "%lf\n", &(solver.dt));
    sscanf(sec->arg_vec_p[2], "%lf\n", &(solver.finaltime));
    if(!procid){
        printf(HEADLINE " cfl: %lf\n", solver.cfl);
        printf(HEADLINE " dt: %lf\n", solver.dt);
        printf(HEADLINE " final time: %lf\n", solver.finaltime);
    }

    /// 4. physical info
    sec = sec_p[4];
    sscanf(sec->arg_vec_p[0], "%lf\n", &(solver.u));
    sscanf(sec->arg_vec_p[1], "%lf\n", &(solver.v));
    sscanf(sec->arg_vec_p[2], "%lf\n", &(solver.viscosity));
    if(!procid){
        printf(HEADLINE " u = %f\n", solver.u);
        printf(HEADLINE " v = %f\n", solver.v);
        printf(HEADLINE " vis = %f\n", solver.viscosity);
    }

    /// 5. LDG parameter
    sec = sec_p[5];
    sscanf(sec->arg_vec_p[0], "%lf\n", &(solver.LDG_parameter[0]));
    sscanf(sec->arg_vec_p[1], "%lf\n", &(solver.LDG_parameter[1]));
    sscanf(sec->arg_vec_p[2], "%lf\n", &(solver.LDG_parameter[2]));
    if(!procid){
        printf(HEADLINE " c11 = %f\n", solver.LDG_parameter[0]);
        printf(HEADLINE " c12 = %f\n", solver.LDG_parameter[1]);
        printf(HEADLINE " c22 = %f\n", solver.LDG_parameter[2]);
    }

    conv_free_section(sec_p);
}

