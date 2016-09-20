#include "SWEDriver2d.h"

/* private function */
void RK45_Coeff(double *, double *, double *);


double SWE_Run2d(PhysDomain2d *phys, SWE_Solver2d *solver, NcFile *outfile){

    /* Runge-Kutta time evaluation coefficient */
    double *rk4a, *rk4b, *rk4c;
    /* allocation */
    rk4a = BuildVector(5);
    rk4b = BuildVector(5);
    rk4c = BuildVector(6);
    RK45_Coeff(rk4a, rk4b, rk4c);

    /* store initial condition */
    int     INTRK, tstep=0;
    int     outstep = 0;
    int     procid = phys->mesh->procid;
    double  time    = 0.0;
    double  ftime   = solver->FinalTime;
    double  dt;  /* delta time */
    double  dtmin = solver->dtmin;
    SWE_StoreVar2d(outfile, phys, outstep++, time);

    double mpitime0 = MPI_Wtime();

    /* time step loop  */
    while (time<ftime){
        dt = SWE_PredictDt2d(phys, solver, 0.3);
        if(dt<dtmin) {dt=dtmin;}
        /* adjust final step to end exactly at FinalTime */
        if (time+dt > ftime) { dt = ftime-time; }

        if(!procid){
            printf("Process:%f, dt:%f\n", time/ftime, dt);
        }

        for (INTRK=1; INTRK<=5; ++INTRK) {
            /* compute rhs of equations */
            const real fdt = (real)dt;
            const real fa = (real)rk4a[INTRK-1];
            const real fb = (real)rk4b[INTRK-1];

            SWE_RHS2d(phys, solver, fa, fb, fdt);
            SLLoc2d(phys, 1.0);
            SWE_PositivePreserving2d(phys, solver);
        }

        time += dt;     /* increment current time */
        tstep++;        /* increment timestep    */
        SWE_StoreVar2d(outfile, phys, outstep++, time);

    }

    double mpitime1 = MPI_Wtime();
    double time_total = mpitime1 - mpitime0;

    return time_total;
}

void RK45_Coeff(double *rk4a, double *rk4b, double *rk4c){
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

