#include <MultiRegions/mr_mesh.h>
#include "conv_driver2d.h"
#include "conv_output.h"
#include "conv_rhs.h"

/* assignment of RK45 parameter */
static void conv_rk_parameter(double *rk4a, double *rk4b, double *rk4c);

/**
 * @brief
 * @details
 * @param [in,out] phys physField object
 */
void conv_run(physField *phys){

    extern conv_solver2d solver;

    double ftime = solver.finaltime;
    double dt = solver.dt;

    double time = 0;
    int    intrk, tstep=0;
    int    counter = 0; // output step
    double rk4a[5], rk4b[5], rk4c[6];

    /* Runge-Kutta time evaluation coefficient */
    conv_rk_parameter(rk4a, rk4b, rk4c);

    /* store initial condition */
    conv_putvar(phys, counter++, time);

    double mpitime0 = MPI_Wtime();
    /* outer time step loop  */
    while (time<ftime){

        /* adjust final step to end exactly at ftime */
        if (time+dt > ftime) { dt = ftime-time; }

        for (intrk=1; intrk<=5; ++intrk) {

            /* compute rhs of equations */
            const real fdt = (real)dt;
            const real fa = (real)rk4a[intrk-1];
            const real fb = (real)rk4b[intrk-1];

            conv_rhs(phys, fa, fb, fdt);
        }

        printf("processing: %f%%\r", time/ftime);
        time += dt;     /* increment current time */
        tstep++;        /* increment timestep    */
        conv_putvar(phys, counter++, time);
    }

    double mpitime1 = MPI_Wtime();
    double elapsetime = mpitime1 - mpitime0;

    MPI_Barrier(MPI_COMM_WORLD);
    const int procid = phys->mesh->procid;
    if(!procid)
        printf("proc: %d,\t time taken: %lg\n", procid, elapsetime);
}

/**
 * @brief set parameters for 4-5 Runge-Kutta method
 */
static void conv_rk_parameter(double *rk4a, double *rk4b, double *rk4c){
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