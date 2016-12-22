//
// Created by li12242 on 12/19/16.
//

#include <MultiRegions/mr_reg.h>
#include "MultiRegions/mr_grid_uniformGrid.h"
#include "LibUtilities/UTest.h"

static int mr_triRegion_test(multiReg *reg, double dt, int verbose);
static int mr_quadRegion_test(multiReg *reg, double dt, int verbose);

int mr_reg_test(int verbose){
    int fail = 0, N = 2;

    extern int Mx, My;

    double clockT1, clockT2;

    /* triangle regions */
    stdCell *shape = sc_create(N, TRIANGLE);
    geoGrid *grid = mr_grid_createUniformGrid_tri(shape, Mx, My, -1, 1, -1, 1, 1);
    clockT1 = MPI_Wtime();
    multiReg *region = mr_reg_create(grid);
    clockT2 = MPI_Wtime();

    fail = mr_triRegion_test(region, clockT2 - clockT1, verbose);

    /* free memory */
    sc_free(shape);
    mr_grid_free(grid);
    mr_reg_free(region);

    /* quadrilateral regions */
    shape = sc_create(N, QUADRIL);
    grid = mr_grid_createUniformGrid_quad(shape, Mx, My, -1, 1, -1, 1);
    clockT1 = MPI_Wtime();
    region = mr_reg_create(grid);
    clockT2 = MPI_Wtime();

    fail = mr_quadRegion_test(region, clockT2 - clockT1, verbose);

    /* free memory */
    sc_free(shape);
    mr_grid_free(grid);
    mr_reg_free(region);

    return fail;
}

static int mr_triRegion_test(multiReg *reg, double dt, int verbose){

    int fail = 0;
    geoGrid *grid = reg->grid;

    const int K = grid->K;
    const int Np = grid->cell->Np;
    const int Nfaces = grid->cell->Nfaces;


    if(verbose){
        /* gen log filename */
        char casename[32] = "mr_triRegion_test";
        FILE *fp = CreateLog(casename, reg->procid, reg->nprocs);

        PrintMatrix2File(fp, "triReg->x", reg->x, K, Np);
        PrintMatrix2File(fp, "triReg->y", reg->y, K, Np);
        PrintMatrix2File(fp, "triReg->drdx", reg->drdx, K, Np);
        PrintMatrix2File(fp, "triReg->drdy", reg->drdy, K, Np);
        PrintMatrix2File(fp, "triReg->dsdx", reg->dsdx, K, Np);
        PrintMatrix2File(fp, "triReg->dsdy", reg->dsdy, K, Np);
        PrintMatrix2File(fp, "triReg->J", reg->J, K, Np);
        PrintMatrix2File(fp, "triReg->nx", reg->nx, K, Nfaces);
        PrintMatrix2File(fp, "triReg->ny", reg->ny, K, Nfaces);
        fclose(fp);
    }

    return fail;
}

static int mr_quadRegion_test(multiReg *reg, double dt, int verbose){

    int fail = 0;
    geoGrid *grid = reg->grid;

    const int K = grid->K;
    const int Np = grid->cell->Np;
    const int Nfaces = grid->cell->Nfaces;


    if(verbose){
        /* gen log filename */
        char casename[32] = "mr_quadRegion_test";
        FILE *fp = CreateLog(casename, reg->procid, reg->nprocs);

        PrintMatrix2File(fp, "quadReg->x", reg->x, grid->K, reg->cell->Np);
        PrintMatrix2File(fp, "quadReg->y", reg->y, grid->K, reg->cell->Np);
        PrintMatrix2File(fp, "triReg->drdx", reg->drdx, K, Np);
        PrintMatrix2File(fp, "triReg->drdy", reg->drdy, K, Np);
        PrintMatrix2File(fp, "triReg->dsdx", reg->dsdx, K, Np);
        PrintMatrix2File(fp, "triReg->dsdy", reg->dsdy, K, Np);
        PrintMatrix2File(fp, "triReg->J", reg->J, K, Np);
        PrintMatrix2File(fp, "triReg->nx", reg->nx, K, Nfaces);
        PrintMatrix2File(fp, "triReg->ny", reg->ny, K, Nfaces);
        fclose(fp);
    }
    return fail;
}
