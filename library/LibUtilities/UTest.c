/**
 * @file
 * Unit test functions
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "LibUtilities.h"

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
int CreateIntMatrixTest(char *message, int **A, int **ExactA, int Nrows, int Ncols){
    double error=0.0, relativeErr;
    double total=0.0;
    int success=0, i, j;
    double **errorMatrix = BuildMatrix(Nrows, Ncols);
    printf("------------Creating %s Test------------\n", message);
    for(i=0;i<Nrows;i++){
        for(j=0;j<Ncols;j++){
            errorMatrix[i][j] = A[i][j] - ExactA[i][j];
            error += fabs( errorMatrix[i][j] );
            total += fabs( ExactA[i][j] );
        }
    }
    relativeErr = error/total;
    printf("Total    Err = %f\n", error);
    printf("Total        = %f\n", total);
    printf("Relative Err = %e\n", relativeErr);
    if(error > TOTALERR | relativeErr > RELATIVEERROR) {
        success = 1; // error flag
        printf("fatal error in %s\n", message);
        PrintMatrix("The error Matrix", errorMatrix, Nrows, Ncols);
    }else{
        printf("test in %s success!\n", message);
    }
    printf("-------------Finish %s Test-------------\n", message);
    DestroyMatrix(errorMatrix);
    return success;
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
int CreateMatrixTest(char *message, double **A, double **ExactA, int Nrows, int Ncols){
    double error=0.0, relativeErr;
    double total=0.0;
    int success=0, i, j;
    double **errorMatrix = BuildMatrix(Nrows, Ncols);
    printf("------------Creating %s Test------------\n", message);
    for(i=0;i<Nrows;i++){
        for(j=0;j<Ncols;j++){
            errorMatrix[i][j] = A[i][j] - ExactA[i][j];
            error += fabs( errorMatrix[i][j] );
            total += fabs( ExactA[i][j] );
        }
    }
    relativeErr = error/total;
    printf("Total    Err = %f\n", error);
    printf("Total        = %f\n", total);
    printf("Relative Err = %e\n", relativeErr);
    if(error > TOTALERR | relativeErr > RELATIVEERROR) {
        success = 1; // error flag
        printf("fatal error in %s\n", message);
        PrintMatrix("The error Matrix", errorMatrix, Nrows, Ncols);
    }else{
        printf("test in %s success!\n", message);
    }
    printf("-------------Finish %s Test-------------\n", message);
    DestroyMatrix(errorMatrix);
    return success;
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
int CreateVectorTest(char *message, double *A, double *ExactA, int Ncols){
    double error=0.0, relativeErr;
    double total=0.0;
    int success=0, i;
    double *errorVector = BuildVector(Ncols);
    printf("------------Creating %s Test------------\n", message);
    for(i=0;i<Ncols;i++){
        errorVector[i] = A[i] - ExactA[i];
        error += fabs( errorVector[i] );
        total += fabs( ExactA[i] );

    }
    relativeErr = error/total;
    printf("Total    Err = %f\n", error);
    printf("Total        = %f\n", total);
    printf("Relative Err = %e\n", relativeErr);
    if(error > TOTALERR | relativeErr > RELATIVEERROR) {
        success = 1; // error flag
        printf("fatal error in %s\n", message);
        PrintVector("The input Vector =", A, Ncols);
        PrintVector("The exact Vector =", ExactA, Ncols);
    }else{
        printf("test in %s success!\n", message);
    }
    printf("-------------Finish %s Test-------------\n", message);

    DestroyVector(errorVector);
    return success;
}

void PrintVector(char *message, double *A, int Ncols){
    int m;
    printf("%s\n", message);
    for(m=0;m<Ncols;++m){
        printf(" %e ", A[m]);
    }
    printf(" \n");
}

void PrintMatrix(char *message, double **A, int Nrows, int Ncols){
    int n,m;

    printf("%s\n", message);
    for(n=0;n<Nrows;++n){
        for(m=0;m<Ncols;++m){
            printf(" %e ", A[n][m]);
        }
        printf(" \n");
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
 * @param [char*] funname   Function name
 * @param [int]   rank      Index of local process
 * @param [int]   nprocs    Number of process
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


