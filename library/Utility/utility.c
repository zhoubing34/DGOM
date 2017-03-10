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

/* sort numbers from small to large */
static int cmp_increase(const void *a, const void *b){
    return (* (int *)a) - (* (int *)b);
}

/**
 * @brief count the number of uique elements in an array
 * @param len length of the array
 * @param list array of integer
 * @return n number of unique elements
 */
int unique_int(int len, int *list){
    if(len ==0) { return 0; }
    qsort(list, (size_t)len, sizeof(int), cmp_increase); // sort the list
    int n=1, i;
    for(i=0;i<(len-1);i++){
        if(list[i+1] != list[i]) { n++; }
    }
    return n;
}