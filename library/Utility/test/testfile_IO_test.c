//
// Created by li12242 on 17/4/14.
//

#include "testfile_IO_test.h"

int testfile_IO_test(int verbose){

    char filename[] = "utility/testfile_IO_test.txt";
    int fail=0,N=10,Nfld=4;
    int **var = matrix_int_create(Nfld, N);
    read_int_file(filename, N, Nfld, var[0], var[1], var[2], var[3]);

    int n,fld,ext[N];
    for(fld=0;fld<Nfld;fld++){
        for(n=0;n<N;n++){ ext[n] = fld+1; }
        fail = vector_int_test(__FUNCTION__, var[fld], ext, N);
    }

    if(verbose){
        int procid, nprocs;
        MPI_Comm_rank(MPI_COMM_WORLD, &procid);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        FILE *fp = create_log(__FUNCTION__, procid, nprocs);
        print_int_matrix2file(fp, "var", var, Nfld, N);
        fclose(fp);
    }

    matrix_int_free(var);
    return fail;
}