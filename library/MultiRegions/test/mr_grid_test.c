//
// Created by li12242 on 12/18/16.
//

#include <StandCell/sc_stdcell.h>
#include <MultiRegions/mr_grid.h>
#include "MultiRegions/mr_grid_create.h"
#include "Utility/UTest.h"

static int mr_triDepartition_test(geoGrid *, double, int verbose);
static int mr_quadDepartition_test(geoGrid *, double, int verbose);

/* number of elements on each coordinate */

int mr_grid_test(int verbose){

    int fail=0, N=2;
    extern int Mx, My;

    double clockT1, clockT2;

    /* triangular geometry grid */
    stdCell *tri = sc_create(N, TRIANGLE);
    clockT1 = MPI_Wtime();
    geoGrid *triGrid = mr_grid_create_uniform_tri(tri, Mx, My, -1, 1, -1, 1, 1);
    clockT2 = MPI_Wtime();
    fail = mr_triDepartition_test(triGrid, clockT2-clockT1, verbose);

    mr_grid_free(triGrid);
    sc_free(tri);

    /* quadrilateral geometry grid */
    stdCell *quad = sc_create(N, QUADRIL);
    clockT1 = MPI_Wtime();
    geoGrid *quadGrid = mr_grid_create_uniform_quad(quad, Mx, My, -1, 1, -1, 1);
    clockT2 = MPI_Wtime();

    fail = mr_quadDepartition_test(quadGrid, clockT2-clockT1, verbose);

    mr_grid_free(quadGrid);
    sc_free(quad);

    return fail;
}


static int mr_triDepartition_test(geoGrid *grid, double dt, int verbose){

    int fail=0, p, Ktotal=0;
    extern int Mx, My;

    /* find number of elements on all processors */
    int Kprocs[grid->nprocs];
    int Klocal = grid->K;
    MPI_Allgather(&Klocal, 1, MPI_INT, Kprocs, 1, MPI_INT, MPI_COMM_WORLD);
    for(p=0;p<grid->nprocs;p++){
        Ktotal += Kprocs[p];
    }
    int total = Mx*My*2;
    fail = IntVector_test("mr_grid_test triangle", &Ktotal, &total, 1, dt);

    if(verbose){
        /* gen log filename */
        char casename[20] = "mr_triGrid_test";
        FILE *fp = CreateLog(casename, grid->procid, grid->nprocs);
        PrintIntMatrix2File(fp, "EToV", grid->EToV, grid->K, grid->cell->Nv);
        fclose(fp);
    }
    return fail;
}

static int mr_quadDepartition_test(geoGrid *grid, double dt, int verbose){

    int fail=0, p, Ktotal=0;
    extern int Mx, My;

    /* find number of elements on all processors */
    int Kprocs[grid->nprocs];
    int Klocal = grid->K;
    MPI_Allgather(&Klocal, 1, MPI_INT, Kprocs, 1, MPI_INT, MPI_COMM_WORLD);
    for(p=0;p<grid->nprocs;p++){
        Ktotal += Kprocs[p];
    }
    int total = Mx*My;
    fail = IntVector_test("mr_grid_test quadrilateral", &Ktotal, &total, 1, dt);

    if(verbose){
        /* gen log filename */
        char casename[20] = "mr_quadGrid_test";
        FILE *fp = CreateLog(casename, grid->procid, grid->nprocs);
        PrintIntMatrix2File(fp, "EToV", grid->EToV, grid->K, grid->cell->Nv);
        fclose(fp);
    }

    return fail;
}
