//
// Created by li12242 on 16/12/27.
//

#include <PhysField/phys_physField.h>
#include "conv_driver2d.h"

void conv_intilization(physField *phys){

    extern conv_solver2d solver;
    double t;
    double dt = 1e6;

    const int K = phys->grid->K;
    const int Np = phys->cell->Np;

    /* initial position */
    const double sigma = 125*1e3/(33*33);
    const double w = 5*M_PI/6;
    const double xc = 0.0;
    const double yc = 0.6;

    int k,n;

    double **x = phys->region->x;
    double **y = phys->region->y;

    /* initial scalar field U = (u, v, c) */
    int sk = 0;
    for(k=0;k<K;++k){
        for(n=0;n<Np;++n){
            const double xt = x[k][n];
            const double yt = y[k][n];
            t = -sigma * ( ( xt - xc )*( xt - xc ) + ( yt - yc )*( yt - yc ) );
            phys->f_Q[sk++] = (real) exp(t); // c field
            phys->f_Q[sk++] = (real)(-w * yt); // flow rate at x-coordinate
            phys->f_Q[sk++] = (real)( w * xt); // flow rate at y-coordinate
        }
    }

    /* time step */
    const double cfl = solver.cfl;
    const int N = solver.N;
    multiReg *reg = phys->region;

    sk = 0;
    for(k=0;k<K;++k){
        double r = reg->len[k];
        for(n=0;n<Np;n++){
            const real u = phys->f_Q[sk++];
            const real v = phys->f_Q[sk++];
            sk++; // jump c field
            double spe = sqrt(u*u+v*v);
            dt = min(dt, r/spe/(N+1));
        }
    }

    double gdt;
    MPI_Allreduce(&dt, &gdt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt = cfl*gdt;

    /* assignment to global variable */
    solver.dt = dt;

    return;
}