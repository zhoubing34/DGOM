//
// Created by li12242 on 17/1/16.
//

#include "arg_section.h"

#define ARG_LEN 120

/**
 * @brief
 * create a pointer as an array of pointers to store the arguments.
 * @param[in] N number of arguments.
 * @return
 * a pointer to pointers of strings.
 */
static char** arg_create(int N){
    char **arg = (char**) calloc((size_t)N, sizeof(char*));
    arg[0] = (char*) calloc( (size_t)N*ARG_LEN, sizeof(char) );
    int i;
    for(i=1;i<N;++i){
        arg[i] = arg[i-1]+ ARG_LEN;
    }
    return arg;
}

/**
 * @brief free the memory of arguments
 */
static void arg_free(char **arg){
    free(arg[0]);
    free(arg);
}


void section_print(arg_section *section_p){
    int i;
    printf("%s", section_p->info_str);
    for(i=0;i<section_p->arg_num;i++){
        printf("%s", section_p->arg_str[i]);
    }
    return;
}

/**
 * @brief create argument section for input parameters.
 *
 * @param[in] info_str information string of argument.
 * @param[in] arg_num number of input argument.
 * @return arg_section argument section structure.
 * @note
 * the `info_str` string should contain `\n` as the last character.
 */
arg_section* section_create(char *info_str, int arg_num){
    arg_section* section_p = (arg_section*) calloc(1, sizeof(arg_section));
    size_t info_lenth = strlen(info_str);
    section_p->info_str = (char*) calloc(info_lenth, sizeof(char));
    strcpy(section_p->info_str, info_str);

    section_p->info_linenum = 0;
    register int i;
    for(i=0;i<info_lenth;i++){
        if( info_str[i] == '\n' )
            section_p->info_linenum++;
    }
    if(info_str[info_lenth-1] != '\n'){
        fprintf(stderr, "section_create (%s): %d\n"
                "error in the info_str - '%s'\n"
                "the last character should be a '\\n'\n",
                __FILE__, __LINE__, info_str);
        exit(-1);
    }
    section_p->arg_num = arg_num;
    section_p->arg_str = arg_create(arg_num);

    return section_p;
}

/**
 * @brief
 * free the memory of argument section structure
 * @param[in] section_p argument section structure
 */
void section_free(arg_section* section_p){
    free(section_p->info_str);
    arg_free(section_p->arg_str);
    free(section_p);
    return;
}

/**
 * @brief write the argument section to file
 * @param section_p argument section structure
 * @param fp handle of file
 */
void section_write_file(arg_section *section_p, FILE *fp){
    int i;
    fprintf(fp, "%s", section_p->info_str);
    for(i=0;i<section_p->arg_num;i++){
        fprintf(fp, "%s\n", section_p->arg_str[i]);
    }
}
/**
 * @brief read argument from file and store it in the structure `section_p`
 * @param section_p argument section strucutre
 * @param fp handle of file
 */
void section_read_file(arg_section *section_p, FILE *fp){
    char buffer[ARG_LEN];
    int i;
    for(i=0;i<section_p->info_linenum;i++){ fgets(buffer, ARG_LEN, fp); }
    for(i=0;i<section_p->arg_num;i++){
        fgets(buffer, ARG_LEN, fp);
        strcpy(section_p->arg_str[i], buffer);
    }
}