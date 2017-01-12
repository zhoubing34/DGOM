/**
 * @file
 * Unit test functions
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "UTest.h"

#define _TOTAL_ERR 1.0e-6
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
 * @brief
 * Compare Matrix variable with exact solution
 *
 * @param [in]  message Test name
 * @param [in]  A
 * @param [in]  ExactA
 * @param [in]  Nrows
 * @param [in]  Ncols
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * fail  | int      | 0 - success; 1 - fail
 */
int IntMatrix_test(char *message,
                   int **A, int **ExactA,
                   int Nrows, int Ncols,
                   double elapsedTime) {

    double error=0.0, relativeErr;
    double total=0.0;
    int fail=0, i, j;

    double **errorMatrix = matrix_double_create(Nrows, Ncols);

    printf(HEADLINE "1 test from %s\n", message);
    for(i=0;i<Nrows;i++){
        for(j=0;j<Ncols;j++){
            errorMatrix[i][j] = A[i][j] - ExactA[i][j];
            error += fabs( errorMatrix[i][j] );
            total += fabs( ExactA[i][j] );
        }
    }
    relativeErr = error/total;

    if(error > _TOTAL_ERR | relativeErr > _RELATIVE_EER) {
        fail = 1; // error flag
        printf(HEADFAIL "1 test failed from %s\n", message);
        printf("Total Err    = %f\n", error);
        printf("Total Value  = %f\n", total);
        printf("Relative Err = %e\n", relativeErr);
        PrintIntMatrix_test("The input Matrix", A, Nrows, Ncols);
        PrintIntMatrix_test("The exact Matrix", ExactA, Nrows, Ncols);
    }else{
        printf(HEADPASS "1 test passed from %s (%f sec)\n", message, elapsedTime);
    }
    matrix_double_free(errorMatrix);
    return fail;
}


/**
 * @brief
 * Compare Matrix variable with exact solution
 *
 * @param [char*]   message
 * @param [double]  A
 * @param [double]  ExactA
 * @param [int]
 * @param [int]
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * fail  | int      | 0 - success; 1 - fail
 */
int Matrix_test(char *message,
                double **A, double **ExactA,
                int Nrows, int Ncols,
                double elapsedTime){

    double error=0.0, relativeErr;
    double total=0.0;
    int fail=0, i, j;

    double **errorMatrix = matrix_double_create(Nrows, Ncols);

    printf(HEADLINE "1 test from %s\n", message);

    for(i=0;i<Nrows;i++){
        for(j=0;j<Ncols;j++){
            errorMatrix[i][j] = A[i][j] - ExactA[i][j];
            error = fabs( errorMatrix[i][j] );
            if (error > _TOTAL_ERR){
                fail = 1; // error flag
                printf(HEADFAIL "1 test failed from %s\n", message);
                printf("element [%d, %d] = %f is different from the exact value %f\n",
                       i+1, j+1, A[i][j], ExactA[i][j]);
                return fail;
            }
            total += fabs( ExactA[i][j] );
        }
    }
    relativeErr = error/total;

    if(relativeErr > _RELATIVE_EER) {
        fail = 1; // error flag
        printf(HEADFAIL "1 test failed from %s\n", message);
        printf("Total Err    = %f\n", error);
        printf("Total Value  = %f\n", total);
        printf("Relative Err = %e\n", relativeErr);

        PrintMatrix_test("The input Matrix", A, Nrows, Ncols);
        PrintMatrix_test("The exact Matrix", ExactA, Nrows, Ncols);
    }else{
        printf(HEADPASS "1 test passed from %s (%f sec)\n", message, elapsedTime);
    }
    matrix_double_free(errorMatrix);
    return fail;
}

/**
 * @brief
 * Compare Matrix variable with exact solution
 *
 * @param [char*]   message
 * @param [double]  A
 * @param [double]  ExactA
 * @param [int]
 * @param [int]
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * success  | int      | 0 - success; 1 - fail
 */
int Vector_test(char *message,
                double *A, double *ExactA, int Ncols,
                double elapsedTime){

    double error=0.0, relativeErr;
    double total=0.0;
    int fail=0, i;

    double *errorVector = vector_double_create(Ncols);

    printf(HEADLINE "1 test from %s\n", message);

    for(i=0;i<Ncols;i++){
        errorVector[i] = A[i] - ExactA[i];
        error += fabs( errorVector[i] );
        total += fabs( ExactA[i] );

    }
    relativeErr = error/total;

    if(error > _TOTAL_ERR | relativeErr > _RELATIVE_EER) {
        fail = 1; // error flag
        printf(HEADFAIL "1 test failed from %s\n", message);

        printf("Total    Err = %f\n", error);
        printf("Total        = %f\n", total);
        printf("Relative Err = %e\n", relativeErr);
        PrintVector_test("The input Vector =", A, Ncols);
        PrintVector_test("The exact Vector =", ExactA, Ncols);
    }else{
        printf(HEADPASS "1 test passed from %s (%f sec)\n", message, elapsedTime);
    }

    vector_double_free(errorVector);
    return fail;
}

int IntVector_test(char *message,
                    int *A, int *ExactA,
                    int Ncols,
                    double elapsedTime){

    double error=0.0, relativeErr;
    double total=0.0;
    int fail=0, i;

    double *errorVector = vector_double_create(Ncols);

    printf(HEADLINE "1 test from %s\n", message);

    for(i=0;i<Ncols;i++){
        errorVector[i] = A[i] - ExactA[i];
        error += fabs( errorVector[i] );
        total += fabs( ExactA[i] );

    }
    relativeErr = error/total;

    if(error > _TOTAL_ERR | relativeErr > _RELATIVE_EER) {
        fail = 1; // error flag
        printf(HEADFAIL "1 test failed from %s\n", message);

        printf("Total    Err = %f\n", error);
        printf("Total        = %f\n", total);
        printf("Relative Err = %e\n", relativeErr);
        PrintIntVector_test("The input Vector =", A, Ncols);
        PrintIntVector_test("The exact Vector =", ExactA, Ncols);
    }else{
        printf(HEADPASS "1 test passed from %s (%f sec)\n", message, elapsedTime);
    }

    vector_double_free(errorVector);
    return fail;
}

/* @brief marco for print vector */
#define _PRINT_VECTOR(message, Ncol, fmt, A) do{\
    int _dim1;\
    printf("%s\n", message);\
    for(_dim1=0;_dim1<Ncol;_dim1++)\
        printf(fmt, A[_dim1]);\
    printf("\n");\
}while(0)

void PrintVector_test(char *message, double *A, int Ncols){
    _PRINT_VECTOR(message, Ncols, " %e ", A);
}
void PrintIntVector_test(char *message, int *A, int Ncols){
    _PRINT_VECTOR(message, Ncols, " %d ", A);
}

/* @brief macros for print matrix */
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


/* @brief macros for write matrix to file */
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
    _WRITE_MATRIX(fp, message, row, col, " %d, ", Mat);
}
void PrintMatrix2File(FILE *fp, char *message, double **Mat, int row, int col){
    _WRITE_MATRIX(fp, message, row, col, " %f, ", Mat);
}

/* @brief marco for write vector to file */
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
    _WRITE_VECTOR(fp, message, len, " %f, ", Mat);
}

/**
 * @brief
 * Create log file to check function routine
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 * @param [in] funname   Function name
 * @param [in] rank      Index of local process
 * @param [in] nprocs    Number of process
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * fig      | FILE*    | file handle
 *
 * @attention
 * Remember to close file with file handle
 *
 */
FILE* CreateLog(char *funname, int rank, int nprocs){

#ifndef DSET_NAME_LEN
#define DSET_NAME_LEN 1024
#endif

    int ret;
    char filename[DSET_NAME_LEN];

    ret = snprintf(filename, DSET_NAME_LEN, "%s%d-%d.txt", funname, rank, nprocs);
    if (ret >= DSET_NAME_LEN) {
        fprintf(stderr, "name too long \n");
        exit(-1);
    }

#undef DSET_NAME_LEN
    FILE *fig = fopen(filename, "w+");

    return fig;
}