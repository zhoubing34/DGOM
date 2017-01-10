//
// Created by li12242 on 16/12/27.
//

#include <PhysField/phys_physField.h>
#include "conv_driver2d.h"
#include "PhysField/phys_add_LDG_solver.h"

static double conv_advectDiff(physField *phys);
static double conv_rotation(physField *phys);

void conv_intilization(physField *phys){

    extern conv_solver2d solver;
    const double cfl = solver.cfl;
    double dt;

    switch (solver.caseid){
        case conv_rotational_convection:
            dt = conv_rotation(phys); break;
        case conv_advection_diffusion:
            dt = conv_advectDiff(phys); break;
        default:
            exit(-1);
    }

    double gdt;
    MPI_Allreduce(&dt, &gdt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt = cfl*gdt;

    /* assignment to global variable */
    solver.dt = dt;
    return;
}

static double conv_advectDiff(physField *phys){

    extern conv_solver2d solver;
    register int k,n;

    const int K = phys->grid->K;
    const int Np = phys->cell->Np;

    double **x = phys->region->x;
    double **y = phys->region->y;

    phys_add_LDG_solver(phys);

    const double xc = -0.5;
    const double yc = -0.5;
    const double miu = solver.viscosity;

    int sk = 0;
    for(k=0;k<K;++k){
        for(n=0;n<Np;++n){
            const double xt = x[k][n];
            const double yt = y[k][n];
            double t = -( ( xt - xc )*( xt - xc ) + ( yt - yc )*( yt - yc ) )/miu;
            phys->f_Q[sk++] = (real) exp(t); // c field
            phys->f_Q[sk++] = (real) solver.u; // flow rate at x-coordinate
            phys->f_Q[sk++] = (real) solver.v; // flow rate at y-coordinate
        }
    }

    sk = 0;
    for(k=0;k<K;++k){
        for(n=0;n<Np;++n){
            phys->viscosity->vis_Q[sk++] = (real) miu; // viscosity parameter for c field
            phys->viscosity->vis_Q[sk++] = (real) 0;
            phys->viscosity->vis_Q[sk++] = (real) 0;
        }
    }


    /* time step */
    double dt = 1e6; // initial time step
    const int N = solver.N;
    double *len = phys->region->len;

    sk = 0;
    for(k=0;k<K;++k){
        double r = len[k]/(N+1);
        for(n=0;n<Np;n++){
            sk++; // jump c field
            const real u = phys->f_Q[sk++];
            const real v = phys->f_Q[sk++];
            double spe = sqrt(u*u+v*v);
            dt = min(dt, r/spe);
            dt = min(dt, r*r/sqrt(miu));
        }
    }

    return dt;
}

static double conv_rotation(physField *phys){

    register int k,n;

    const int K = phys->grid->K;
    const int Np = phys->cell->Np;

    double **x = phys->region->x;
    double **y = phys->region->y;

    /* initial position */
    const double sigma = 125*1e3/(33*33);
    const double w = 5*M_PI/6;
    const double xc = 0.0;
    const double yc = 0.6;

    /* initial scalar field U = (u, v, c) */
    int sk = 0;
    for(k=0;k<K;++k){
        for(n=0;n<Np;++n){
            const double xt = x[k][n];
            const double yt = y[k][n];
            double t = -sigma * ( ( xt - xc )*( xt - xc ) + ( yt - yc )*( yt - yc ) );
            phys->f_Q[sk++] = (real) exp(t); // c field
            phys->f_Q[sk++] = (real)(-w * yt); // flow rate at x-coordinate
            phys->f_Q[sk++] = (real)( w * xt); // flow rate at y-coordinate
        }
    }

    extern conv_solver2d solver;
    /* time step */
    double dt = 1e6; // initial time step
    const int N = solver.N;
    double *len = phys->region->len;

    sk = 0;
    for(k=0;k<K;++k){
        double r = len[k]/(N+1);
        for(n=0;n<Np;n++){
            sk++; // jump c field
            const real u = phys->f_Q[sk++];
            const real v = phys->f_Q[sk++];
            double spe = sqrt(u*u+v*v);
            dt = min(dt, r/spe);
        }
    }

    return dt;
}