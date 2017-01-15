//
// Created by li12242 on 17/1/13.
//

#ifndef DGOM_SWE_READ_COMMAND_H
#define DGOM_SWE_READ_COMMAND_H

typedef enum {
    swe_command_help,
    swe_command_create_input,
    swe_command_run,
}swe_command_type;

swe_command_type swe_read_command(int argc, char **argv);

void swe_print_help();

#endif //DGOM_SWE_READ_COMMAND_H
