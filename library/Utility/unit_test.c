#include "unit_test.h"

#define _TOTAL_ERR 1.0e-10

/**
 * @brief get parameters from commandline
 * @param [in]     argc
 * @param [in]     argv
 * @param [in,out] isverbose
 * @param [in,out] ishelp
 *
 */
void test_command(int argc, char **argv, int *ishelp, int *isverbose){
    *ishelp = 0;
    *isverbose = 0;

    if(argc>1){
        int i;
        for(i=0; i<argc; i++){
            if(!(memcmp(argv[i], "-help", 5)) )
                *ishelp = 1;
            else if(!(memcmp(argv[i], "-verbose", 8)) )
                *isverbose = 1;
        }
    }
    return;
}

#define check_matrix_err(fmt) \
do{ \
    if (abserr > _TOTAL_ERR){                                                   \
        fail = 1;                                                               \
        printf(HEADFAIL "1 test failed from %s\n", message);                    \
        printf("element [%d, %d] = ", i, j);                                    \
        printf(fmt, A[i][j]);                                                   \
        printf(" is different from the exact value ");                          \
        printf(fmt, A_ext[i][j]);                                               \
        printf("\n");                                                           \
        return fail;                                                            \
    }                                                                           \
}while(0)

/**
 * @brief compare matrix with the exact solution
 * @param message
 * @param A
 * @param A_ext
 * @param Nrows
 * @param Ncols
 * @return
 */
int matrix_int_test(const char *message, int **A, int **A_ext, int Nrows, int Ncols)
{
    int fail=0,i,j;
    double tmp, abserr;
    for(i=0;i<Nrows;i++){
        for(j=0;j<Ncols;j++){
            tmp = A[i][j] - A_ext[i][j];
            abserr = abs( tmp );
            check_matrix_err("%d");
        }
    }
    return fail;
}

/**
 * @brief compare matrix variable with the exact solution
 * @param message
 * @param A
 * @param A_ext
 * @param Nrows
 * @param Ncols
 * @return
 */
int matrix_double_test(const char *message, double **A, double **A_ext, int Nrows, int Ncols)
{
    int fail=0,i,j;
    double tmp, abserr;
    for(i=0;i<Nrows;i++){
        for(j=0;j<Ncols;j++){
            tmp = A[i][j] - A_ext[i][j];
            abserr = abs( tmp );
            check_matrix_err("%f");
        }
    }
    return fail;
}

#define check_vector_err(fmt) \
do{ \
    if (abserr > _TOTAL_ERR){                                                   \
        fail = 1;                                                               \
        printf(HEADFAIL "1 test failed from %s\n", message);                    \
        printf("element [%d] = ", i);                                           \
        printf(fmt, A[i]);                                                      \
        printf(" is different from the exact value ");                          \
        printf(fmt, A_ext[i]);                                                  \
        printf("\n");                                                           \
        return fail;                                                            \
    }                                                                           \
}while(0)

/**
 * @brief check vector with the exact solution
 * @param message
 * @param A
 * @param A_ext
 * @param Ncols
 * @return
 */
int vector_double_test(const char *message, double *A, double *A_ext, int Ncols)
{

    double abserr,tmp;
    int fail=0,i;

    for(i=0;i<Ncols;i++){
        tmp = A[i] - A_ext[i];
        abserr = abs( tmp );
        check_vector_err("%f");
    }
    return fail;
}
/**
 * @brief check vector with the exact solution
 * @param message
 * @param A
 * @param A_ext
 * @param Ncols
 * @return
 */
int vector_int_test(const char *message, int *A, int *A_ext, int Ncols)
{
    double abserr,tmp;
    int fail=0,i;

    for(i=0;i<Ncols;i++){
        tmp = A[i] - A_ext[i];
        abserr = abs( tmp );
        check_vector_err("%d");
    }
    return fail;
}

/**
 * @brief marco to print vector value
 */
#define _PRINT_VECTOR(message, Ncol, fmt, A) \
do{                                         \
    int _dim1;                              \
    printf("%s\n", message);                \
    for(_dim1=0;_dim1<Ncol;_dim1++)         \
        printf(fmt, A[_dim1]);              \
    printf("\n");                           \
}while(0)

/**
 * @brief print vector of double type
 * @param message message heading
 * @param A vector
 * @param Ncols vector length
 */
void print_double_vector(char *message, double *A, int Ncols){
    _PRINT_VECTOR(message, Ncols, " %e ", A);
}
void print_int_vector(char *message, int *A, int Ncols){
    _PRINT_VECTOR(message, Ncols, " %d ", A);
}

/**
 * @brief macros to print matrix
 * */
#define _PRINT_MATRIX(message, Nrow, Ncol, fmt, A)  \
do{                                                 \
    int _dim1, _dim2;                               \
    printf("%s\n", message);                        \
    for(_dim1=0;_dim1<Nrow;_dim1++){                \
        for(_dim2=0;_dim2<Ncol;_dim2++)             \
            printf(fmt, A[_dim1][_dim2]);           \
        printf("\n");                               \
    }                                               \
}while(0)

void print_double_matrix(char *message, double **A, int Nrows, int Ncols){
    _PRINT_MATRIX(message, Nrows, Ncols, " %e ", A);
}
void print_int_matrix(char *message, int **A, int Nrows, int Ncols){
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


void print_int_matrix2file(FILE *fp, char *message, int **Mat, int row, int col){
    _WRITE_MATRIX(fp, message, row, col, " %d ", Mat);
}
void print_double_matrix2file(FILE *fp, char *message, double **Mat, int row, int col){
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

void print_int_vector2file(FILE *fp, char *message, int *Mat, int len){
    _WRITE_VECTOR(fp, message, len, " %d, ", Mat);
}
void print_double_vector2file(FILE *fp, char *message, double *Mat, int len){
    _WRITE_VECTOR(fp, message, len, " %f ", Mat);
}

/**
 * @brief create log file to check function routine
 *
 * @param[in] name name string
 * @param[in] rank index of local process
 * @param[in] nprocs number of process
 * @return fig file handle of log file.
 * @author li12242, Tianjin University, li12242@tju.edu.cn
 */
FILE* create_log(const char *name, int rank, int nprocs){
    int ret;
    char filename[MAX_NAME_LENGTH];
    ret = snprintf(filename, MAX_NAME_LENGTH, "%s.%d-%d.txt", name, rank, nprocs);
    if (ret >= MAX_NAME_LENGTH) {
        fprintf(stderr, "%s (%s): %d\nthe function name %s is too long \n",
                __FUNCTION__, __FILE__, __LINE__, name);
        exit(-1);
    }
    FILE *fig = fopen(filename, "ws+");
    return fig;
}