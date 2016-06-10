#include "Convection2d/Convection2d.h"

void ConvectionRun2d(Mesh *mesh, Ncfile * outfile, double FinalTime, double dt){
    double time = 0;
    int    INTRK, tstep=0;
    int    counter = 0;

    PutVar(outfile, counter++, time, mesh);

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

#ifdef LIMIT
            DisDetector(mesh); /* discontinuity detector */
            LimiterBJ2d(mesh);
#endif
        }

        time += dt;     /* increment current time */
        tstep++;        /* increment timestep    */

    }

    double mpitime1 = MPI_Wtime();
    double time_total = mpitime1 - mpitime0;

    PutVar(outfile, counter++, time, mesh);

    double flopsV = p_Np*p_Np*12 + p_Np*13; /* V3 */
    double flopsS = p_Nfp*p_Nfaces*21 + p_Np*(p_Nfaces*p_Nfp*6 + 3);
    double flopsR = p_Np*p_Nfields*4;
    int Kloc = mesh->K;

    MPI_Barrier(MPI_COMM_WORLD);
    printf("proc: %d,\t order: %d,\t time taken: %lg,\t MNUPS: %lg,\t GFLOPS: %lg\n",
           mesh->procid, p_N, time_total, p_N*mesh->K*5*(tstep-1)/time_total*1e-6,
           5*(tstep-1)*( (flopsV+flopsS+flopsR)*((double)Kloc/(1.e9*time_total))));
}
