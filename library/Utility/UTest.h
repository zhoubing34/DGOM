#ifndef DGOM_UTEST_H
#define DGOM_UTEST_H

#include "LibUtilities.h"
#include <sys/time.h>

#define HEADSTART   "\n\033[32m[   RUNING   ]\033[0m " // module reminder
#define HEADLINE 	"\033[32m[------------]\033[0m " // test reminder
#define HEADPASS  	"\033[32m[   PASSED   ]\033[0m " // test pass indicator
#define HEADFAIL 	"\033[31m[   FAILED   ]\033[0m " // test fail indicator
#define HEADEND     "\033[32m[============]\033[0m " // test end indicator

/* UTest.c */
void UTest_Command(int argc, char** argv, int *ishelp, int *isverbose);

void PrintVector_test(char *message, double *A, int Ncols);
void PrintIntVector_test(char *message, int *A, int Ncols);
void PrintMatrix_test(char *message, double **A, int Nrows, int Ncols);
void PrintIntMatrix_test(char *message, int **A, int Nrows, int Ncols);
void PrintIntMatrix2File(FILE *fp, char *message, int **Mat, int row, int col);
void PrintMatrix2File(FILE *fp, char *message, double **Mat, int row, int col);
void PrintIntVector2File(FILE *fp, char *message, int *Mat, int len);
void PrintVector2File(FILE *fp, char *message, double *Mat, int len);

int IntMatrix_test(char *message, int **A, int **ExactA, int Nrows, int Ncols, double elapsedTime);
int Matrix_test(char *message, double **A, double **ExactA, int Nrows, int Ncols, double elapsedTime);
int Vector_test(char *message, double *A, double *ExactA, int Ncols, double elapsedTime);
int IntVector_test(char *message, int *A, int *ExactA, int Ncols, double elapsedTime);

FILE* CreateLog(char *funname, int nprocs, int rank);

#endif //DGOM_UTEST_H
