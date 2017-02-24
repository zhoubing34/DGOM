//
// Created by li12242 on 17/2/24.
//
#include "utility.h"

/**
 * @brief
 * Transfer string to integer.
 * @details
 * Usages:
 *      str2int(argv[1], &N , "Wrong degree input");
 *
 * @param[in] str the input string
 * @param[out] N  output integer
 * @param[in] errmessage error message
 */
void str2int(char *str, int *N, char* errmessage){
    if (sscanf(str,"%d",N) !=1) {
        fprintf(stderr, "%s:%s \n", errmessage, str);
        exit(-1);
    }
}

/**
 * @brief
 * transfer string to double.
 * @param[in] str character string
 * @param[out] scal output double variable
 * @param[in] errmessage error message
 */
void str2double(char *str, double *scal, char* errmessage){
    if (sscanf(str,"%lf", scal)!=1) {
        fprintf(stderr, "%s:%s \n", errmessage, str);
        exit(-1);
    }
}
