/**
 * @file
 * unit test functions
 * @author li12242, Tianjin University, li12242@tju.edu.cn
 */
#ifndef DGOM_UTEST_H
#define DGOM_UTEST_H

#include "utility.h"

#define HEADSTART   "\n\033[32m[   RUNING   ]\033[0m " ///< module reminder
#define HEADLINE 	"\033[32m[------------]\033[0m " ///< test reminder
#define HEADPASS  	"\033[32m[   PASSED   ]\033[0m " ///< test pass indicator
#define HEADFAIL 	"\033[31m[   FAILED   ]\033[0m " ///< test fail indicator
#define HEADEND     "\033[32m[============]\033[0m " ///< test end indicator

void test_command(int argc, char **argv, int *ishelp, int *isverbose);

void print_double_vector(char *message, double *A, int Ncols);
void print_int_vector(char *message, int *A, int Ncols);
void print_double_matrix(char *message, double **A, int Nrows, int Ncols);
void print_int_matrix(char *message, int **A, int Nrows, int Ncols);
void print_int_matrix2file(FILE *fp, char *message, int **Mat, int row, int col);
void print_double_matrix2file(FILE *fp, char *message, double **Mat, int row, int col);
void print_int_vector2file(FILE *fp, char *message, int *Mat, int len);
void print_double_vector2file(FILE *fp, char *message, double *Mat, int len);

int matrix_int_test(const char *message, int **A, int **A_ext, int Nrows, int Ncols);
int matrix_double_test(const char *message, double **A, double **A_ext, int Nrows, int Ncols);
int vector_double_test(const char *message, double *A, double *A_ext, int Ncols);
int vector_int_test(const char *message, int *A, int *A_ext, int Ncols);

FILE* create_log(const char *name, int nprocs, int rank);

#endif //DGOM_UTEST_H
