#include "Convection2d/Convection2d.h"
#include <mpi.h>

void ConvectionRun2d(Mesh *mesh, double FinalTime, double dt){
    double time = 0;
    int    INTRK, tstep=0;

    double mpitime0 = MPI_Wtime();

    /* outer time step loop  */
    while (time<FinalTime){

        /* adjust final step to end exactly at FinalTime */
        if (time+dt > FinalTime) { dt = FinalTime-time; }

        for (INTRK=1; INTRK<=5; ++INTRK) {

            /* compute rhs of equations */
            const float fdt = dt;
            const float fa = (float)mesh->rk4a[INTRK-1];
            const float fb = (float)mesh->rk4b[INTRK-1];

            ConvectionRHS2d(mesh, fa, fb, fdt);
        }

        time += dt;     /* increment current time */
        tstep++;        /* increment timestep    */

        Write2TestFile(mesh, time);
    }

    double mpitime1 = MPI_Wtime();

    double time_total = mpitime1-mpitime0;

    MPI_Barrier(MPI_COMM_WORLD);
    printf("proc: %d,\t order: %d,\t time taken: %lg\n", mesh->procid, p_N, time_total);
    
}
