#include <MultiRegions/mr_mesh.h>
#include "conv_driver2d.h"
#include "PhysField/phys_physField.h"
#include "conv_output.h"
#include "conv_rhs.h"

/* global variable for RK45 parameter */
double rk4a[5];
double rk4b[5];
double rk4c[5];

/* assignment of RK45 parameter */
void conv_rk_parameter();

/**
 * @brief
 * @details
 * @param [in,out] phys physField object
 */
void conv_run(physField *phys){

    extern conv_solver2d solver;

    double ftime = solver.finalTime;
    double dt = solver.dt;

    double time = 0;
    int    intrk, tstep=0;
    int    counter = 0; // output step
    extern double rk4a[5], rk4b[5], rk4c[5];

    /* Runge-Kutta time evaluation coefficient */
    conv_rk_parameter();

    /* store initial condition */
    conv_putvar(phys, counter++, time);

    double mpitime0 = MPI_Wtime();
    /* outer time step loop  */
    while (time<ftime){

        /* adjust final step to end exactly at ftime */
        if (time+dt > ftime) { dt = ftime-time; }

        for (intrk=1; intrk<=5; ++intrk) {

            /* compute rhs of equations */
            const float fdt = (float)dt;
            const float fa = (float)rk4a[intrk-1];
            const float fb = (float)rk4b[intrk-1];

            conv_rhs(phys, fa, fb, fdt);

        }

        time += dt;     /* increment current time */
        tstep++;        /* increment timestep    */
        conv_putvar(phys, counter++, time);

    }

    double mpitime1 = MPI_Wtime();
    double elapsetime = mpitime1 - mpitime0;

    int Np      = phys->cell->Np;
    int Nfaces  = phys->cell->Nfaces;
    int Nfp     = phys->cell->Nfp;
    int Nfields = phys->Nfield;
    int N       = phys->cell->N;

    double flopsV = Np*Np*12 + Np*13; /* V3 */
    double flopsS = Nfp*Nfaces*21 + Np*(Nfaces*Nfp*6 + 3);
    double flopsR = Np*Nfields*4;
    int Kloc = phys->grid->K;
    int procid = phys->mesh->procid;
    int nprocs = phys->mesh->nprocs;

    MPI_Barrier(MPI_COMM_WORLD);
    printf("proc: %d,\t order: %d,\t time taken: %lg,\t MNUPS: %lg,\t GFLOPS: %lg\n",
           procid, N, elapsetime, N*Kloc*5*(tstep-1)/elapsetime*1e-6,
           5*(tstep-1)*( (flopsV+flopsS+flopsR)*((double)Kloc/(1.e9*elapsetime))));
}

/**
 * @brief set parameters for 4-5 Runge-Kutta method
 */
void conv_rk_parameter(){
    extern double rk4a[5], rk4b[5], rk4c[5];
    /* low storage RK coefficients */
    rk4a[0] =              0.0;
    rk4a[1] =  -567301805773.0 / 1357537059087.0;
    rk4a[2] = -2404267990393.0 / 2016746695238.0;
    rk4a[3] = -3550918686646.0 / 2091501179385.0;
    rk4a[4] = -1275806237668.0 /  842570457699.0;
    rk4b[0] =  1432997174477.0 /  9575080441755.0;
    rk4b[1] =  5161836677717.0 / 13612068292357.0;
    rk4b[2] =  1720146321549.0 /  2090206949498.0;
    rk4b[3] =  3134564353537.0 /  4481467310338.0;
    rk4b[4] =  2277821191437.0 / 14882151754819.0;
    rk4c[0] =              0.0;
    rk4c[1] =  1432997174477.0 / 9575080441755.0;
    rk4c[2] =  2526269341429.0 / 6820363962896.0;
    rk4c[3] =  2006345519317.0 / 3224310063776.0;
    rk4c[4] =  2802321613138.0 / 2924317926251.0;
    rk4c[5] =              1.0;
}