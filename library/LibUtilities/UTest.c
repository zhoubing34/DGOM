/**
 * @file
 * Unit test functions
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "UTest.h"

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

    double **errorMatrix = Matrix_create(Nrows, Ncols);

    printf(HEADLINE "1 test from %s\n", message);
    for(i=0;i<Nrows;i++){
        for(j=0;j<Ncols;j++){
            errorMatrix[i][j] = A[i][j] - ExactA[i][j];
            error += fabs( errorMatrix[i][j] );
            total += fabs( ExactA[i][j] );
        }
    }
    relativeErr = error/total;

    if(error > TOTALERR | relativeErr > RELATIVEERROR) {
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
    Matrix_free(errorMatrix);
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

    double **errorMatrix = Matrix_create(Nrows, Ncols);

    printf(HEADLINE "1 test from %s\n", message);

    for(i=0;i<Nrows;i++){
        for(j=0;j<Ncols;j++){
            errorMatrix[i][j] = A[i][j] - ExactA[i][j];
            error += fabs( errorMatrix[i][j] );
            total += fabs( ExactA[i][j] );
        }
    }
    relativeErr = error/total;

    if(error > TOTALERR | relativeErr > RELATIVEERROR) {
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
    Matrix_free(errorMatrix);
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
                double *A, double *ExactA,
                int Ncols,
                double elapsedTime){

    double error=0.0, relativeErr;
    double total=0.0;
    int fail=0, i;

    double *errorVector = Vector_create(Ncols);

    printf(HEADLINE "1 test from %s\n", message);

    for(i=0;i<Ncols;i++){
        errorVector[i] = A[i] - ExactA[i];
        error += fabs( errorVector[i] );
        total += fabs( ExactA[i] );

    }
    relativeErr = error/total;

    if(error > TOTALERR | relativeErr > RELATIVEERROR) {
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

    Vector_free(errorVector);
    return fail;
}

int IntVector_test(char *message,
                    int *A, int *ExactA,
                    int Ncols,
                    double elapsedTime){

    double error=0.0, relativeErr;
    double total=0.0;
    int fail=0, i;

    double *errorVector = Vector_create(Ncols);

    printf(HEADLINE "1 test from %s\n", message);

    for(i=0;i<Ncols;i++){
        errorVector[i] = A[i] - ExactA[i];
        error += fabs( errorVector[i] );
        total += fabs( ExactA[i] );

    }
    relativeErr = error/total;

    if(error > TOTALERR | relativeErr > RELATIVEERROR) {
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

    Vector_free(errorVector);
    return fail;
}

void PrintVector_test(char *message, double *A, int Ncols){
    int m;
    printf("%s\n", message);
    for(m=0;m<Ncols;++m){
        printf(" %e ", A[m]);
    }
    printf(" \n");
}

void PrintIntVector_test(char *message, int *A, int Ncols){
    int m;
    printf("%s\n", message);
    for(m=0;m<Ncols;++m){
        printf(" %d ", A[m]);
    }
    printf(" \n");
}

void PrintMatrix_test(char *message, double **A, int Nrows, int Ncols){
    int n,m;

    printf("%s\n", message);
    for(n=0;n<Nrows;++n){
        for(m=0;m<Ncols;++m){
            printf(" %e ", A[n][m]);
        }
        printf(" \n");
    }
}

void PrintIntMatrix_test(char *message, int **A, int Nrows, int Ncols){
    int n,m;

    printf("%s\n", message);
    for(n=0;n<Nrows;++n){
        for(m=0;m<Ncols;++m){
            printf(" %d ", A[n][m]);
        }
        printf(" \n");
    }
}

void PrintIntMatrix2File(FILE *fp, char *message, int **Mat, int row, int col){
    fprintf(fp, "%s = \n", message);
    int n,m;
    for(n=0;n<row;++n){
        for(m=0;m<col;++m){
            fprintf(fp, " %d, ", Mat[n][m]);
        }
        fprintf(fp, " \n");
    }
}

void PrintIntVector2File(FILE *fp, char *message, int *Mat, int len){
    fprintf(fp, "%s = \n", message);
    int n;
    for(n=0;n<len;++n){
        fprintf(fp, " %d, ", Mat[n]);
    }
    fprintf(fp, "\n");
}

void PrintVector2File(FILE *fp, char *message, double *Mat, int len){
    fprintf(fp, "%s = \n", message);
    int n;
    for(n=0;n<len;++n){
        fprintf(fp, " %f, ", Mat[n]);
    }
    fprintf(fp, "\n");
}

void PrintMatrix2File(FILE *fp, char *message, double **Mat, int row, int col){
    fprintf(fp, "%s = \n", message);
    int n,m;
    for(n=0;n<row;++n){
        for(m=0;m<col;++m){
            fprintf(fp, " %f, ", Mat[n][m]);
        }
        fprintf(fp, " \n");
    }
}

void SaveMatrix(char *filename, double **A, int Nrows, int Ncols){
    int n,m;
    FILE *fp = fopen(filename, "w");
    for(n=0;n<Nrows;++n){
        for(m=0;m<Ncols;++m){
            fprintf(fp, " %g ", A[n][m]);
        }
        fprintf(fp, " \n");
    }

    fclose(fp);
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

    FILE *fig = fopen(filename, "w");

    return fig;
}


