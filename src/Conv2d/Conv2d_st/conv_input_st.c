//
// Created by li12242 on 17/3/28.
//

#include "conv2d_st.h"

/* read command line */
static Conv_Run_Type conv_read_command(int argc, char **argv);

//==============================function for conv2d_st====================================//
/**
 * @brief create argument section structure for standard convection problem.
 * @return argument section structure.
 */
static arg_section** conv_arg_section_create_st(){
    arg_section** section_p = (arg_section**) calloc(SEC_NUM, sizeof(arg_section*));

    int ind = 0;
    char probleminfo[] = HEADEND "DGOM: 2d convection problem\n"
    HEADLINE "case indicator (1 parameter)\n"
    HEADLINE "    1. case indicator  |-- 0. rotational convection\n"
    HEADLINE "                       |-- 1. advection-diffusion\n"
    HEADLINE "                       |-- 2. pure advection\n";
    int var_num = 1;
    section_p[ind++] = section_create(probleminfo, var_num);

    char scinfo[] = HEADEND "standard element info (2 parameters)\n"
    HEADLINE "    1. element type  |--- 2. triangle\n"
    HEADLINE "                     |--- 3. quadrilateral\n"
    HEADLINE "    2. order of polynomial;\n";
    var_num = 2;
    section_p[ind++] = section_create(scinfo, var_num);

    char meshinfo[] = HEADEND "mesh info (3 parameters)\n"
    HEADLINE "    1. num of elements in x direction;\n"
    HEADLINE "    2. num of elements in y direction;\n";
    var_num = 2;
    section_p[ind++] = section_create(meshinfo, var_num);

    char timeinfo[] = HEADEND "time info (3 parameters)\n"
    HEADLINE "    1. CFL number;\n"
    HEADLINE "    2. dt;\n"
    HEADLINE "    3. final time;\n";
    var_num = 3;
    section_p[ind++] = section_create(timeinfo, var_num);

    char physinfo[] = HEADEND "parameters for advection-diffusion case\n"
    HEADLINE "    1. flow rate u for x direction;\n"
    HEADLINE "    2. flow rate v for x direction;\n"
    HEADLINE "    3. viscosity miu;\n";
    var_num = 3;
    section_p[ind++] = section_create(physinfo, var_num);

    char outputinfo[] = HEADEND "output file\n"
    HEADLINE "    1. output file name;\n"
    HEADLINE "    2. output time interval;\n";
    var_num = 2;
    section_p[ind] = section_create(outputinfo, var_num);

    return section_p;
}

/** @brief create input file. */
static void conv_create_inputfile_st(char *filename){
    arg_section **sec_p = conv_arg_section_create_st();
    int n;
    for(n=0;n<SEC_NUM;n++) {section_print(sec_p[n]);}

    FILE *fp = fopen(filename, "w");
    for(n=0;n<SEC_NUM;n++){ section_write_file(fp, sec_p[n]); }
    fclose(fp);
    conv_arg_section_free(sec_p);
    return;
}

/**
 * @brief read parameters from input command or file
 * @param [in] argc number of conmand
 * @param [in] argv pointers to the command strings
 */
void conv_input_st(int argc, char **argv){
    /* read from command */
    Conv_Run_Type run_type = conv_read_command(argc, argv);

    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    /* const strings */
    static char helpinfo[] = HEADEND "DGOM:\n" HEADLINE "2d convection problem\n"
            HEADLINE "Optional features:\n"
            HEADLINE "   -help          print help information;\n"
            HEADLINE "   -create_input  create input file;\n"
            HEADLINE "   -run           run with input file;\n"
            HEADLINE "Example usages:\n"
            HEADLINE "   mpirun -n 2 -host localhost ./conv2d_st -create_input conv_parameter.inc\n"
            HEADLINE "\n"
            HEADLINE "   mpirun -n 2 -host localhost ./conv2d_st -run conv_parameter.inc\n";

    if(run_type == CONV_HELP){
        if(!procid) { printf("%s", helpinfo); } exit(0);
    }else if(run_type == CONV_CREATE_INPUT){
        if( (!procid) & (argc < 2) ) { fprintf(stderr, "Unknown input filename\n"); exit(-1); }
        if(!procid) { conv_create_inputfile_st(argv[2]); } exit(0);
    }else if(run_type == CONV_RUN){
        if( (!procid) & (argc < 2) ) { fprintf(stderr, "Unknown input filename\n"); exit(-1); }
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
    arg_section **sec_p = conv_arg_section_create_st();
    /* open file and read section */
    FILE *fp = fopen(filename, "r");
    int n;
    for(n=0;n<SEC_NUM;n++){
        section_read_file(fp, sec_p[n]);
    }
    fclose(fp);
    return sec_p;
}

/**
 * @brief read all the command line arguments and return the Conv_Run_Type.
 * @param argc number of command argument;
 * @param argv pointers to command argument;
 * @return
 * run type.
 */
static Conv_Run_Type conv_read_command(int argc, char **argv){
    Conv_Run_Type command_type = CONV_HELP;

    char help_str[] = "-help";
    char pre_str[]  = "-create_input";
    char run_str[]  = "-run";

    register int i;
    for(i=0;i<argc;i++){
        if(!(memcmp(argv[i], help_str, strlen(help_str))) ){
            command_type = CONV_HELP;
        }else if(!(memcmp(argv[i], pre_str, strlen(pre_str) ))){
            command_type = CONV_CREATE_INPUT;
        }else if(!(memcmp(argv[i], run_str, strlen(run_str) ))){
            command_type = CONV_RUN;
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