//
// Created by li12242 on 17/1/13.
//

#include <string.h>
#include "swe_read_command.h"
#include "Utility/unit_test.h"

swe_command_type swe_read_command(int argc, char **argv){

    swe_command_type command_type = swe_command_run;

    char help_str[] = "-help";
    char create_str[] = "-preprocess";

    int i;
    for(i=0;i<argc;i++){
        if(!(memcmp(argv[i], help_str, strlen(help_str))) ){
            command_type = swe_command_help; break;
        }else if(!(memcmp(argv[i], create_str, strlen(create_str))) ){
            command_type = swe_command_create_input; break;
        }
    }

    return command_type;
}


void swe_print_help(){

    char helpinfo[] = HEADEND "DGOM:\n" HEADLINE "2d swe solver\n"
    HEADLINE "Optional features:\n"
    HEADLINE "   -help         print help information\n"
    HEADLINE "   -preprocess   generate input file\n"
    HEADLINE "Example usages:\n"
    HEADLINE "   mpirun -n 2 -host localhost ./swe2d -preprocess\n"
    HEADEND "\n";

    printf("%s", helpinfo);
}