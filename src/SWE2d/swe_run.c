#include <PhysField/pf_phys.h>
#include "swe_driver2d.h"
#include "PhysField/pf_cellMean.h"
#include "swe_output.h"
#include "PhysField/Limiter/pf_limit_BJ2d.h"
#include "swe_rhs.h"
#include "PhysField/pf_openbc.h"

#define DEBUG 0

/* private function */
void swe_rk_parameter(double *rk4a, double *rk4b, double *rk4c);
double swe_time_interval(swe_solver *solver);
void swe_ppreserve(swe_solver *solver);

static void swe_h2eta(swe_solver *solver){
    physField *phys = solver->phys;
    const int K = phys->grid->K;
    const int Np = phys->cell->Np;
    const int Nfield = phys->Nfield;
    dg_real *f_Q = phys->f_Q;

    register int k,n,sk;
    for(k=0;k<K;k++){
        for(n=0;n<Np;n++){
            sk = k*Np+n;
            dg_real bot = solver->bot[sk];
            f_Q[sk*Nfield] -= bot;
        }
    }
}

static void swe_eta2h(swe_solver *solver){
    physField *phys = solver->phys;
    const int K = phys->grid->K;
    const int Np = phys->cell->Np;
    const int Nfield = phys->Nfield;
    dg_real *f_Q = phys->f_Q;

    register int k,n,sk;
    for(k=0;k<K;k++){
        for(n=0;n<Np;n++){
            sk = k*Np+n;
            dg_real bot = solver->bot[sk];
            f_Q[sk*Nfield] += bot;
        }
    }
}

void swe_run(swe_solver *solver){
    /* Runge-Kutta time evaluation coefficient */
    double rk4a[5], rk4b[5], rk4c[6];
    swe_rk_parameter(rk4a, rk4b, rk4c);

    physField *phys = solver->phys;

    /* store loop condition */
    int     INTRK, tstep=0;
    int     procid = phys->mesh->procid;
    double  time    = 0.0;
    double  tloc = 0.0;
    double  ftime   = solver->ftime;
    double  dt;  /* delta time */
    /* save initial condition */
    swe_save_var(solver, tstep++, time);
    pf_limit_BJ2d(phys, 1.0);
    swe_ppreserve(solver);

    double mpitime0 = MPI_Wtime();

    /* time step loop  */
    while (time<ftime){
        /* save result */
        swe_save_var(solver, tstep++, time);

        /* calculate time interval */
        dt = swe_time_interval(solver);

        //if(dt<1e-4) dt = 1e-4;
        /* adjust final step to end exactly at FinalTime */
        if (time+dt > ftime) { dt = ftime-time; }

         if(!procid){
             printf("Process:%f, dt:%f\r", time/ftime, dt);
         }

        for (INTRK=1; INTRK<=5; ++INTRK) {
            /* compute rhs of equations */
            const dg_real fdt = (dg_real)dt;
            const dg_real fa = (dg_real)rk4a[INTRK-1];
            const dg_real fb = (dg_real)rk4b[INTRK-1];
            pf_set_openbc(phys, time+rk4c[INTRK-1]*dt, time_interp_linear);
            swe_rhs(solver, fa, fb, fdt);
            //swe_h2eta(solver);
            pf_limit_BJ2d(phys, 1.0);
            //swe_eta2h(solver);
            //pf_vert_limit(phys);
            swe_ppreserve(solver);
        }
        /* increment current time */
        time += dt;
    }
    double mpitime1 = MPI_Wtime();
    double elapsetime = mpitime1 - mpitime0;

    // last
    swe_save_var(solver, tstep++, time);

    if(!procid)
        printf("proc: %d, time taken: %lg\n", procid, elapsetime);

    return;
}

void swe_rk_parameter(double *rk4a, double *rk4b, double *rk4c){
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

/**
 * @brief
 * Predict the wave speed and give local delta time
 *
 * @details
 * The wave speed is given as \f[ c = \sqrt{gh}+\left| \mathbf{u} \right| \f],
 * and the delta time is derived \f[ dt = dx/dt \f], while \f$ dx \f$ is the
 * characteristic length.
 *
 * @param[PhysDomain2d*] phys
 * @param[SWE_Solver2d*] solver
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * dt   | double | delta time
 *
 */
double swe_time_interval(swe_solver *solver){

    double dt   = 1e4;
    const double gra  = solver->gra;
    const double hcrit = solver->hcrit;

    physField *phys = solver->phys;
    const int Np = phys->cell->Np;
    const int K = phys->grid->K;
    const int Nfield = phys->Nfield;

    register int k,n;
    double gdt;
    dg_real h,u,v;

    for(k=0;k<K;k++){
        double len = phys->region->len[k];
        double c = 1e-10;
        int isdry = 0;
        for(n=0;n<Np;n++){
            int ind = (k*Np+n)*Nfield;
            h = phys->f_Q[ind++];
            u = phys->f_Q[ind++]; u /= h;
            v = phys->f_Q[ind  ]; v /= h;

            if(h<hcrit){
                isdry = 1;
            }
            c = max(c, sqrt(gra*h)+sqrt(u*u+v*v));
        }
        if(isdry){ continue; } // jump this cell
        dt = min(dt, len/(double)c);
#if 0
        int procid;
        MPI_Comm_rank(MPI_COMM_WORLD, &procid);
        if(!procid) printf("k=%d, hmean=%f, c=%f, dt=%f\n",k,hmean,c,dt);
#endif
    }
    /* gather all the delta time in all process */
    dt = dt*solver->cfl;
    MPI_Allreduce(&dt, &gdt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return gdt;
}


/**
 * @brief
 * Positive-preserving operator from Xing and Shu (2010).
 *
 * @details
 * The positive operator will reconstruct the conserve variables in each element
 * to eliminate the negative water depth.
 *
 * @param[in,out] solver SWE_Solver2d pointer
 */
void swe_ppreserve(swe_solver *solver){

    physField *phys = solver->phys;

    const int Nfield = phys->Nfield;
    const int Np = phys->cell->Np;
    const int K = phys->grid->K;
    const double hcrit = solver->hcrit;

    register int i,k,ind;

    pf_cellMean(phys);

    for(k=0;k<K;k++){
        ind = k*Nfield; /* index of k-th cell */
        dg_real hmean  = phys->c_Q[ind];
        dg_real qxmean = phys->c_Q[ind+1];
        dg_real qymean = phys->c_Q[ind+2];
        /* correct negative water depth */
        if(hmean<0.0){
            for(i=0;i<Np;i++) {
                ind = (k*Np + i)*Nfield;
                phys->f_Q[ind  ] -= hmean;
                phys->f_Q[ind+1]  = 0.0;
                phys->f_Q[ind+2]  = 0.0;
            }
            hmean = 0.0;
            qxmean = 0.0;
            qymean = 0.0;
        }

        ind = k*Np*Nfield; /* index of first node */
        dg_real hmin = phys->f_Q[ind];
        /* compute minimum water depth */
        for(i=1;i<Np;i++){
            ind += Nfield;
            hmin = min(hmin, phys->f_Q[ind]);
        }
        /* positive operator */
        dg_real theta;
        if(hmean > hmin){ /* in case for `hmean = hmin` */
            theta = min(1, hmean/(hmean-hmin));
        }else{ /* for hmean = hmin */
            theta = 0.0;
        }
        /* reconstruction */
        for(i=0;i<Np;i++){
            ind = (k*Np + i)*Nfield;
            phys->f_Q[ind  ] = (phys->f_Q[ind  ] - hmean )*theta + hmean;
            phys->f_Q[ind+1] = (phys->f_Q[ind+1] - qxmean)*theta + qxmean;
            phys->f_Q[ind+2] = (phys->f_Q[ind+2] - qymean)*theta + qymean;

            if(phys->f_Q[ind]<hcrit){
                phys->f_Q[ind+1] = 0.0;
                phys->f_Q[ind+2] = 0.0;
            }
        }
    }

    return;
}
