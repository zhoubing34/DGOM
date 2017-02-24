//
// Created by li12242 on 17/1/16.
//
/**
 * @file
 * 参数文件处理
 * @brief
 * pre-process structure for input files
 * @details
 * 1. `section_create` return a pointer to a new arg_section;
 * 2. `section_write_file` write the arg_section to the input file;
 * 3. `section_read_file` read argument from input file to arg_section.
 * 4. `section_free` deallocate the memory of arg_section.
 */

#ifndef DGOM_ARG_SECTION_H
#define DGOM_ARG_SECTION_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct{
    int info_linenum; ///< number of lines in information string
    char *info_str; ///< strings of information
    int arg_num; ///< number of argument
    char **arg_vec_p; ///< string of argument
} arg_section;

void section_print(arg_section *section_p);
arg_section* section_create(char *info_str, int arg_num);
void section_free(arg_section* section_p);
void section_write_file(arg_section *section_p, FILE *fp);
void section_read_file(arg_section *section_p, FILE *fp);

#endif //DGOM_ARG_SECTION_H
