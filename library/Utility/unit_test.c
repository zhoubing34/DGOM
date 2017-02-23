/**
 * @file Unit test functions
 *
 * @author li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "unit_test.h"

#define _TOTAL_ERR 1.0e-10
#define _RELATIVE_EER 1.0e-8

/**
 * @brief Get configures from commandline and return settings
 * @param [in]     argc
 * @param [in]     argv
 * @param [in,out] isverbose
 * @param [in,out] ishelp
 *
 */
void UTest_Command(int argc, char** argv, int *ishelp, int *isverbose){
    int i;
    *ishelp = 0;
    *isverbose = 0;

    if(argc>1){
        for(i=0; i<argc; i++){
            if(!(memcmp(argv[i], "-help", 5)) )
                *ishelp = 1;
            else if(!(memcmp(argv[i], "-verbose", 8)) )
                *isverbose = 1;
        }
    }
    return;
}

/**
 * @brief Compare matrix with exact solution
 * @param message
 * @param A
 * @param ExactA
 * @param Nrows
 * @param Ncols
 * @return
 */
int matrix_int_test(const char *message, int **A, int **ExactA, int Nrows, int Ncols)
{
    int fail=0,i,j;
    double tmp, abserr;
    for(i=0;i<Nrows;i++){
        for(j=0;j<Ncols;j++){
            tmp = A[i][j] - ExactA[i][j];
            abserr = fabs( tmp );
            if (abserr > _TOTAL_ERR){
                fail = 1; // error flag
                printf(HEADFAIL "1 test failed from %s\n", message);
                printf("element [%d, %d] = %d is different from the exact value %d\n",
                       i, j, A[i][j], ExactA[i][j]);
                return fail;
            }
        }
    }
    return fail;
}


/**
 * @brief Compare matrix variable with exact solution
 * @param message
 * @param A
 * @param ExactA
 * @param Nrows
 * @param Ncols
 * @return
 */
int matrix_double_test(const char *message, double **A, double **ExactA, int Nrows, int Ncols)
{
    int fail=0,i,j;
    double tmp, abserr;
    for(i=0;i<Nrows;i++){
        for(j=0;j<Ncols;j++){
            tmp = A[i][j] - ExactA[i][j];
            abserr = fabs( tmp );
            if (abserr > _TOTAL_ERR){
                fail = 1; // error flag
                printf(HEADFAIL "1 test failed from %s\n", message);
                printf("element [%d, %d] = %f is different from the exact value %f\n",
                       i, j, A[i][j], ExactA[i][j]);
                return fail;
            }
        }
    }
    return fail;
}

/**
 * @brief
 * @param message
 * @param A
 * @param ExactA
 * @param Ncols
 * @return
 */
int vector_double_test(const char *message, double *A, double *ExactA, int Ncols)
{

    double abserr,tmp;
    int fail=0,i;

    for(i=0;i<Ncols;i++){
        tmp = A[i] - ExactA[i];
        abserr = fabs( tmp );
        if (abserr > _TOTAL_ERR){
            fail = 1; // error flag
            printf(HEADFAIL "1 test failed from %s\n", message);
            printf("element [%d] = %f is different from the exact value %f\n",
                   i, A[i], ExactA[i]);
            return fail;
        }
    }
    return fail;
}
/**
 * @brief
 * @param message
 * @param A
 * @param ExactA
 * @param Ncols
 * @return
 */
int vector_int_test(const char *message, int *A, int *ExactA, int Ncols)
{
    double abserr,tmp;
    int fail=0,i;

    for(i=0;i<Ncols;i++){
        tmp = A[i] - ExactA[i];
        abserr = fabs( tmp );
        if (abserr > _TOTAL_ERR){
            fail = 1; // error flag
            printf(HEADFAIL "1 test failed from %s\n", message);
            printf("element [%d] = %d is different from the exact value %d\n",
                   i, A[i], ExactA[i]);
            return fail;
        }
    }
    return fail;
}

/**
 * @brief marco to print vector value
 */
#define _PRINT_VECTOR(message, Ncol, fmt, A) do{\
    int _dim1;\
    printf("%s\n", message);\
    for(_dim1=0;_dim1<Ncol;_dim1++)\
        printf(fmt, A[_dim1]);\
    printf("\n");\
}while(0)

/**
 * @brief print vector of double type
 * @param message message heading
 * @param A vector
 * @param Ncols vector length
 */
void PrintVector(char *message, double *A, int Ncols){
    _PRINT_VECTOR(message, Ncols, " %e ", A);
}
void PrintIntVector(char *message, int *A, int Ncols){
    _PRINT_VECTOR(message, Ncols, " %d ", A);
}

/**
 * @brief macros to print matrix
 * */
#define _PRINT_MATRIX(message, Nrow, Ncol, fmt, A) do{\
    int _dim1, _dim2;\
    printf("%s\n", message);\
    for(_dim1=0;_dim1<Nrow;_dim1++){\
        for(_dim2=0;_dim2<Ncol;_dim2++)\
            printf(fmt, A[_dim1][_dim2]);\
        printf("\n");\
    }\
}while(0)

void PrintMatrix_test(char *message, double **A, int Nrows, int Ncols){
    _PRINT_MATRIX(message, Nrows, Ncols, " %e ", A);
}
void PrintIntMatrix_test(char *message, int **A, int Nrows, int Ncols){
    _PRINT_MATRIX(message, Nrows, Ncols, " %d ", A);
}


/**
 * @brief macros for write matrix to file
 * */
#define _WRITE_MATRIX(fp, message, Nrow, Ncol, fmt, A) do{\
    int _dim1, _dim2;\
    fprintf(fp, "%s=\n", message);\
    for(_dim1=0;_dim1<Nrow;_dim1++){\
        for(_dim2=0;_dim2<Ncol;_dim2++)\
            fprintf(fp, fmt, A[_dim1][_dim2]);\
        fprintf(fp, "\n");\
    }\
}while(0)


void PrintIntMatrix2File(FILE *fp, char *message, int **Mat, int row, int col){
    _WRITE_MATRIX(fp, message, row, col, " %d ", Mat);
}
void PrintMatrix2File(FILE *fp, char *message, double **Mat, int row, int col){
    _WRITE_MATRIX(fp, message, row, col, " %f, ", Mat);
}

/**
 * @brief marco for write vector to file
 * */
#define _WRITE_VECTOR(fp, message, Ncol, fmt, A) do{\
    int _dim1;\
    fprintf(fp, "%s=\n", message);\
    for(_dim1=0;_dim1<Ncol;_dim1++)\
        fprintf(fp, fmt, A[_dim1]);\
    fprintf(fp, "\n");\
}while(0)

void PrintIntVector2File(FILE *fp, char *message, int *Mat, int len){
    _WRITE_VECTOR(fp, message, len, " %d, ", Mat);
}
void PrintVector2File(FILE *fp, char *message, double *Mat, int len){
    _WRITE_VECTOR(fp, message, len, " %f ", Mat);
}

/**
 * @brief Create log file to check function routine
 * @author li12242, Tianjin University, li12242@tju.edu.cn
 *
 * @param[in] funname   Function name
 * @param[in] rank      Index of local process
 * @param[in] nprocs    Number of process
 * @return fig file handle of log file.
 */
FILE* create_log(const char *funname, int rank, int nprocs){

#ifndef DSET_NAME_LEN
#define DSET_NAME_LEN 1024
#endif
    int ret;
    char filename[DSET_NAME_LEN];

    ret = snprintf(filename, DSET_NAME_LEN, "%s%d-%d.txt", funname, rank, nprocs);
    if (ret >= DSET_NAME_LEN) {
        fprintf(stderr, "%s (%s): %d\nthe function name %s is too long \n",
                __FUNCTION__, __FILE__, __LINE__, funname);
        exit(-1);
    }

#undef DSET_NAME_LEN
    FILE *fig = fopen(filename, "w+");
    return fig;
}