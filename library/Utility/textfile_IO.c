//
// Created by li12242 on 17/4/14.
//

#include "textfile_IO.h"
#include <stdarg.h>

/**
 * @brief read the integer arrays from input text file.
 * @param filename name of input file;
 * @param Nlen length of each array;
 * @param Nfield number of array;
 * @param ... pointer to each array;
 */
int read_int_file(char *filename, int Nlen, int Nfield, ...){

    FILE *fp;
    if( (fp = fopen(filename, "r")) == NULL ){
        fprintf(stderr, "%s (%d): Unable to open file %s.\n",
                __FUNCTION__,__LINE__,filename);
    }
    /* read file title */
    char buffer[MAX_NAME_LENGTH];
    fgets(buffer, MAX_NAME_LENGTH, fp);
    int n,fld;
    sscanf(buffer, "%d %d", &n, &fld);
    /* check file title */
    if( n != Nlen ){
        fprintf(stderr, "%s (%d): The length of file %s (%d) is incorrect, expected to be %d\n",
                __FUNCTION__, __LINE__, filename, n, Nlen);
        return -1;
    }
    if( fld != Nfield ){
        fprintf(stderr, "%s (%d): The field of file %s (%d) is incorrect, expected to be %d\n",
                __FUNCTION__, __LINE__, filename, fld, Nfield);
        return -1;
    }

    /* initilize Nfield pointer */
    va_list intptr;
    va_start(intptr, Nfield);
    int *var[Nfield];
    for(fld=0;fld<Nfield;fld++){
        var[fld] = va_arg(intptr, int*);
    }
    int temp;
    /* read file */
    for(n=0;n<Nlen;n++){
        fscanf(fp, "%d", &temp); //read index
        for(fld=0;fld<Nfield;fld++){
            fscanf(fp, "%d", var[fld]+n);
        }
    }

    fclose(fp);
    va_end(intptr);
    return 0;
}

/**
 * @brief read the double arrays from input text file.
 * @param filename name of input file;
 * @param Nlen length of each array;
 * @param Nfield number of array;
 * @param ... pointer to each array;
 */
int read_double_file(char *filename, int Nlen, int Nfield, ...){

    FILE *fp;
    if( (fp = fopen(filename, "r")) == NULL ){
        fprintf(stderr, "%s (%d): Unable to open file %s.\n",
                __FUNCTION__,__LINE__,filename);
    }
    /* read file title */
    char buffer[MAX_NAME_LENGTH];
    fgets(buffer, MAX_NAME_LENGTH, fp);
    int n,fld;
    sscanf(buffer, "%d %d", &n, &fld);
    /* check file title */
    if( n != Nlen ){
        fprintf(stderr, "%s (%d): The length of file %s (%d) is incorrect, expected to be %d\n",
                __FUNCTION__, __LINE__, filename, n, Nlen);
        return -1;
    }
    if( fld != Nfield ){
        fprintf(stderr, "%s (%d): The field of file %s (%d) is incorrect, expected to be %d\n",
                __FUNCTION__, __LINE__, filename, fld, Nfield);
        return -1;
    }

    /* initilize Nfield pointer */
    va_list intptr;
    va_start(intptr, Nfield);
    dg_real *var[Nfield];
    for(fld=0;fld<Nfield;fld++){
        var[fld] = va_arg(intptr, dg_real*);
    }
    int temp;
    /* read file */
    for(n=0;n<Nlen;n++){
        fscanf(fp, "%d", &temp); //read index
        for(fld=0;fld<Nfield;fld++){
            fscanf(fp, "%lf", var[fld]+n);
        }
    }

    fclose(fp);
    va_end(intptr);
    return 0;
}