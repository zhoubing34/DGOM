//
// Created by li12242 on 16/12/27.
//

#include "conv_driver2d.h"
#include "PhysField/pf_add_LDG_solver.h"
#include "PhysField/pf_init_file.h"

#define DEBUG 0

static double conv_advect_diff_init(physField *phys);
static double conv_rotation_init(physField *phys);
static double conv_userset_init(physField *phys);

void conv_intilization(physField *phys){

    extern conv_solver2d solver;
    const double cfl = solver.cfl;
    double dt;

    switch (solver.caseid){
        case conv_rotational_convection:
            dt = conv_rotation_init(phys); break;
        case conv_advection_diffusion:
            dt = conv_advect_diff_init(phys); break;
        case conv_userset:
            dt = conv_userset_init(phys); break;
        default:
            fprintf(stderr, "%s: %d\nUnknown test case %d\n",
                    __FUNCTION__, __LINE__, solver.caseid);
            exit(-1);
    }

    double gdt;
    MPI_Allreduce(&dt, &gdt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt = cfl*gdt;

    /* assignment to global variable */
    solver.dt = dt;

    return;
}

static double conv_userset_init(physField *phys){
    extern conv_solver2d solver;
    pf_init_file2d(phys, solver.casename);
    double dt = 1e6; // initial time step
    const int N = solver.N;
    double *len = phys->region->len;

    const int K = phys->grid->K;
    const int Np = phys->cell->Np;
    const double miu = solver.viscosity;

    int k,n, sk = 0;
    for(k=0;k<K;++k){
        double r = len[k]/(N+1);
        for(n=0;n<Np;n++){
            sk++; // jump c field
            const dg_real u = phys->f_Q[sk++];
            const dg_real v = phys->f_Q[sk++];
            double spe = sqrt(u*u+v*v);
            dt = min(dt, r/spe);
            if(miu > EPS){ dt = min(dt, r*r/sqrt(miu)); }
#if DEBUG
            int procid = phys->mesh->procid;
            if(!procid)
                printf("k=%d, n=%d, r=%f, spe=%f, mu=%f, dt=%f\n", k,n,r,spe,miu,dt);
#endif
        }
    }
    return dt;
}

static double conv_advect_diff_init(physField *phys){

    extern conv_solver2d solver;
    register int k,n;

    const int K = phys->grid->K;
    const int Np = phys->cell->Np;

    double **x = phys->region->x;
    double **y = phys->region->y;

    pf_add_LDG_solver(phys);

    const double xc = -0.5;
    const double yc = -0.5;
    const double miu = solver.viscosity;
    const double sigma = 125*1e3/(33*33);

    int sk = 0;
    for(k=0;k<K;++k){
        for(n=0;n<Np;++n){
            const double xt = x[k][n];
            const double yt = y[k][n];
            double t;
            if(miu > DIFF_THRESHOLD){
                t = -( ( xt - xc )*( xt - xc ) + ( yt - yc )*( yt - yc ) )/miu;
            }else{
                t = -( ( xt - xc )*( xt - xc ) + ( yt - yc )*( yt - yc ) )*sigma;
            }

            phys->f_Q[sk++] = (dg_real) exp(t); // c field
            phys->f_Q[sk++] = (dg_real) solver.u; // flow rate at x-coordinate
            phys->f_Q[sk++] = (dg_real) solver.v; // flow rate at y-coordinate
        }
    }

    sk = 0;
    for(k=0;k<K;++k){
        for(n=0;n<Np;++n){
            phys->viscosity->vis_Q[sk++] = (dg_real) miu; // viscosity parameter for c field
            phys->viscosity->vis_Q[sk++] = (dg_real) 0;
            phys->viscosity->vis_Q[sk++] = (dg_real) 0;
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
            const dg_real u = phys->f_Q[sk++];
            const dg_real v = phys->f_Q[sk++];
            double spe = sqrt(u*u+v*v);
            dt = min(dt, r/spe);
            dt = min(dt, r*r/sqrt(miu));
#if DEBUG
            int procid = phys->mesh->procid;
            if(!procid)
                printf("k=%d, n=%d, r=%f, spe=%f, mu=%f, dt=%f\n", k,n,r,spe,miu,dt);
#endif
        }
    }
    return dt;
}

static double conv_rotation_init(physField *phys){

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
            phys->f_Q[sk++] = (dg_real) exp(t); // c field
            phys->f_Q[sk++] = (dg_real)(-w * yt); // flow rate at x-coordinate
            phys->f_Q[sk++] = (dg_real)( w * xt); // flow rate at y-coordinate
        }
    }

    extern conv_solver2d solver;
    solver.viscosity = 0; // eliminate viscosity
    /* time step */
    double dt = 1e6; // initial time step
    const int N = solver.N;
    double *len = phys->region->len;

    sk = 0;
    for(k=0;k<K;++k){
        double r = len[k]/(N+1);
        for(n=0;n<Np;n++){
            sk++; // jump c field
            const dg_real u = phys->f_Q[sk++];
            const dg_real v = phys->f_Q[sk++];
            double spe = sqrt(u*u+v*v);
            dt = min(dt, r/spe);
        }
    }

    return dt;
}