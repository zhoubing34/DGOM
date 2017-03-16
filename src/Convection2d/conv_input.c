//
// Created by li12242 on 16/12/27.
//

#include "conv_driver.h"

/* read input file */
Conv_Case_Type conv_read_type();
/* create input file */
static void conv_create_inputfile();
/* read command line */
static Conv_Run_Type conv_read_command(int argc, char **argv);
/* create argument section */
arg_section** conv_arg_section_create();

/**
 * @brief read parameters from input command or file
 * @param [in] argc number of conmand
 * @param [in] argv pointers to the command strings
 */
Conv_Case_Type conv_input(int argc, char **argv){
    /* read from command */
    Conv_Run_Type run_type = conv_read_command(argc, argv);

    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    /* const strings */
    char helpinfo[] = HEADEND "DGOM:\n" HEADLINE "2d convection problem\n"
            HEADLINE "Optional features:\n"
            HEADLINE "   -help         print help information\n"
            HEADLINE "   -preprocess   generate input file\n"
            HEADLINE "Example usages:\n"
            HEADLINE "   mpirun -n 2 -host localhost ./convection2d -preprocess\n";

    Conv_Case_Type case_type = conv_userset; // default
    if(run_type == conv_help){
        if(!procid) { printf("%s", helpinfo); } exit(0);
    }else if(run_type == conv_preprocss){
        if(!procid) { conv_create_inputfile(); } exit(0);
    }else if(run_type == conv_execute){
        case_type = conv_read_type();
    }
    return case_type;
}
/**
 * @brief read command line and get the run type.
 * @param argc number of command argument.
 * @param argv pointers to command argument.
 * @return
 * run type.
 */
static Conv_Run_Type conv_read_command(int argc, char **argv){
    Conv_Run_Type command_type = conv_execute;

    char help_str[] = "-help";
    char pre_str[] = "-preprocess";

    register int i;
    for(i=0;i<argc;i++){
        if(!(memcmp(argv[i], help_str, strlen(help_str))) ){
            command_type = conv_help;
        }else if(!(memcmp(argv[i], pre_str, strlen(pre_str) ))){
            command_type = conv_preprocss;
        }
    }
    return command_type;
}

/**
 * @brief create input file
 */
static void conv_create_inputfile(){
    arg_section **sec_p = conv_arg_section_create();
    int n;
    for(n=0;n<SEC_NUM;n++) {section_print(sec_p[n]);}

    FILE *fp = fopen(INPUT_FILE_NAME, "w");
    for(n=0;n<SEC_NUM;n++){
        section_write_file(sec_p[n], fp);
    }
    fclose(fp);
    conv_arg_section_free(sec_p);
    return;
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
/**
 * @brief create argument section structure for convection problem.
 * @return argument section structure.
 */
static arg_section** conv_arg_section_create(){
    arg_section** section_p = (arg_section**) calloc(SEC_NUM, sizeof(arg_section*));

    int ind = 0;
    char probleminfo[] = HEADEND "DGOM: 2d convection problem\n"
            HEADLINE "case indicator (1 parameter)\n"
            HEADLINE "    1. case indicator  |-- 1. rotational convection\n"
            HEADLINE "                       |-- 2. advection-diffusion\n"
            HEADLINE "                       |-- 3. user specific\n";
    int var_num = 1;
    section_p[ind++] = section_create(probleminfo, var_num);

    char scinfo[] = HEADEND "standard element info (2 parameters)\n"
            HEADLINE "    1. element type  |--- 2. triangle\n"
            HEADLINE "                     |--- 3. quadrilateral\n"
            HEADLINE "    2. order of polynomial\n";
    var_num = 2;
    section_p[ind++] = section_create(scinfo, var_num);

    char meshinfo[] = HEADEND "mesh info (3 parameters)\n"
            HEADLINE "    1. num of elements in x direction\n"
            HEADLINE "    2. num of elements in y direction\n"
            HEADLINE "    3. case name\n";
    var_num = 3;
    section_p[ind++] = section_create(meshinfo, var_num);

    char timeinfo[] = HEADEND "time info (3 parameters)\n"
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
 * @brief read argument section and return case type.
 */
Conv_Case_Type conv_read_type(){
    arg_section **sec_p = conv_arg_section_create();
    // read section from file
    FILE *fp = fopen(INPUT_FILE_NAME, "r");
    section_read_file(sec_p[0], fp);
    fclose(fp);

    Conv_Case_Type case_type;

    /// 0. section: case info
    arg_section *sec = sec_p[0];
    sscanf(sec->arg_vec_p[0], "%d\n", &(case_type));
    return case_type;
}

arg_section** conv_read_input(){
    arg_section **sec_p = conv_arg_section_create();
    // read section from file
    FILE *fp = fopen(INPUT_FILE_NAME, "r");
    int n;
    for(n=0;n<SEC_NUM;n++){
        section_read_file(sec_p[n], fp);
    }
    fclose(fp);
    return sec_p;
}

