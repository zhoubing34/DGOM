//
// Created by li12242 on 17/3/19.
//

#include "nc_library_test.h"
#include "Utility/nc_library.h"

int nc_file_read_from_file_test(int verbose){
    int procid, nprocs, fail = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    char filename[] = "SWE2d/TsuamiRunup2d/TsuamiRunup.obc4.nc";
    NC_File *file = nc_file_read_from_file(filename, procid, nprocs);

    if(!procid) {
        nc_file_print(file);
        nc_var_print(file->var_vec_p[2]);
    }
    nc_file_free(file);
    return fail;
}
