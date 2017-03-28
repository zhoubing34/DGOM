//
// Created by li12242 on 16/12/27.
//

#include "conv2d.h"


/* read command line */
static Conv_Run_Type conv_read_command(int argc, char **argv);

/* const strings */
static char helpinfo[] = HEADEND "DGOM:\n" HEADLINE "2d convection problem\n"
        HEADLINE "Optional features:\n"
        HEADLINE "   -help          print help information;\n"
        HEADLINE "   -create_input  create input file;\n"
        HEADLINE "   -run           run with input file;\n"
        HEADLINE "Example usages:\n"
        HEADLINE "   mpirun -n 2 -host localhost ./conv2d -create_input conv_parameter.inc\n"
        HEADLINE "\n"
        HEADLINE "   mpirun -n 2 -host localhost ./conv2d -run conv_parameter.inc\n";

//==============================function for conv2d====================================//

/**
 * @brief create argument section structure for user specific convection problem.
 * @return argument section structure.
 */
static arg_section** conv_arg_section_create(){
    arg_section** section_p = (arg_section**) calloc(SEC_NUM, sizeof(arg_section*));

    /// 0. case info
    int ind = 0;
    char probleminfo[] = HEADEND "DGOM: 2d convection problem\n"
            HEADLINE "case name (1 parameter)\n"
            HEADLINE "    1. case indicator  |-- 1. casename\n";
    int var_num = 1;
    section_p[ind++] = section_create(probleminfo, var_num);

    /// 1. cell info
    char scinfo[] = HEADEND "standard element info (2 parameters)\n"
            HEADLINE "    1. element type  |--- 2. triangle\n"
            HEADLINE "                     |--- 3. quadrilateral\n"
            HEADLINE "    2. order of polynomial;\n";
    var_num = 2;
    section_p[ind++] = section_create(scinfo, var_num);

    /// 2. obc info
    char meshinfo[] = HEADEND "mesh info (3 parameters)\n"
            HEADLINE "    1. open boundary condition file;\n";
    var_num = 1;
    section_p[ind++] = section_create(meshinfo, var_num);

    /// 3. time info
    char timeinfo[] = HEADEND "time info (3 parameters)\n"
            HEADLINE "    1. CFL number;\n"
            HEADLINE "    2. dt;\n"
            HEADLINE "    3. final time;\n";
    var_num = 3;
    section_p[ind++] = section_create(timeinfo, var_num);

    /// 4. initial condition
    char physinfo[] = HEADEND "physical parameter for advection-diffusion case (3 parameters)\n"
            HEADLINE "    1. initial condition file (C and u,v);\n";
    var_num = 1;
    section_p[ind++] = section_create(physinfo, var_num);

    /// 5. viscosity info
    char LDGinfo[] = HEADEND "LDG parameter for viscosity term (3 parameters)\n"
            HEADLINE "    1. C11 (0.0 recommend);\n"
            HEADLINE "    2. C12 (0.5 recommend);\n"
            HEADLINE "    3. C22;\n";
    var_num = 3;
    section_p[ind] = section_create(LDGinfo, var_num);

    return section_p;
}

/** @brief create input file. */
static void conv_create_inputfile(char *filename){
    arg_section **sec_p = conv_arg_section_create();
    int n;
    for(n=0;n<SEC_NUM;n++) {section_print(sec_p[n]);}

    FILE *fp = fopen(filename, "w");
    for(n=0;n<SEC_NUM;n++){ section_write_file(fp, sec_p[n]); }
    fclose(fp);
    conv_arg_section_free(sec_p);
    return;
}


/**
 * @brief read parameters from input command or file.
 * @param [in] argc number of conmand;
 * @param [in] argv pointers to the command strings;
 */
void conv_input(int argc, char **argv){
    /* read from command */
    Conv_Run_Type run_type = conv_read_command(argc, argv);

    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    if(run_type == conv_help){
        if(!procid) { printf("%s", helpinfo); } exit(0);
    }else if(run_type == conv_create_input){
        if(!procid) { conv_create_inputfile(argv[2]); } exit(0);
    }else if(run_type == conv_execute){
        extern Conv_Solver solver;
        strcpy(solver.filename, argv[2]); // set the parameter file
        printf(HEADLINE " input file: %s\n", solver.filename);
    }
    return;
}

/**
 * @brief read parameters from input file.
 * @return
 */
arg_section** conv_read_inputfile(char *filename){
    arg_section **sec_p = conv_arg_section_create();
    /* open file and read section */
    FILE *fp = fopen(filename, "r");
    int n;
    for(n=0;n<SEC_NUM;n++){
        section_read_file(fp, sec_p[n]);
    }
    fclose(fp);
    return sec_p;
}

//==============================utility functions====================================//
/**
 * @brief read all the command line arguments and return the Conv_Run_Type.
 * @param argc number of command argument;
 * @param argv pointers to command argument;
 * @return
 * run type.
 */
static Conv_Run_Type conv_read_command(int argc, char **argv){
    Conv_Run_Type command_type = conv_execute; // default

    char help_str[] = "-help";
    char pre_str[] = "-create_input";

    register int i;
    for(i=0;i<argc;i++){
        if(!(memcmp(argv[i], help_str, strlen(help_str))) ){
            command_type = conv_help;
        }else if(!(memcmp(argv[i], pre_str, strlen(pre_str) ))){
            command_type = conv_create_input;
        }
    }
    return command_type;
}

/**
 * @brief free argument section memory.
 * @param section_p argument section strucutre.
 */
void conv_arg_section_free(arg_section **section_p){
    int n;
    for(n=0;n<SEC_NUM;n++)
        section_free(section_p[n]);
    free(section_p);
    return;
}

