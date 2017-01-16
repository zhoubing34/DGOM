//
// Created by li12242 on 17/1/16.
//

#ifndef DGOM_ARG_SECTION_H
#define DGOM_ARG_SECTION_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct{
    int info_linenum; ///< number of lines
    char *info_str; ///< strings of info
    int arg_num; ///< argument number
    char **arg_str; ///< argument str
} arg_section;

void section_print(arg_section *section_p);
arg_section* section_create(char *info_str, int arg_num);
void section_free(arg_section* section_p);
void section_write_file(arg_section *section_p, FILE *fp);
void section_read_file(arg_section *section_p, FILE *fp);

#endif //DGOM_ARG_SECTION_H
