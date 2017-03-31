#include "swe_lib.h"

#define DEBUG 0

/* private function */
static void swe_rk_parameter(double *rk4a, double *rk4b, double *rk4c);
static double swe_time_interval(dg_phys *phys);
static void swe_ppreserve(dg_phys *phys);

void swe_run(){

    /* Runge-Kutta time evaluation coefficient */
    double rk4a[5], rk4b[5], rk4c[6];
    swe_rk_parameter(rk4a, rk4b, rk4c);

    extern SWE_Solver solver;
    dg_phys *phys = solver.phys;
    /* store loop condition */
    int     INTRK, tstep=0;
    int     procid = dg_mesh_procid(dg_phys_mesh(phys));
    double  time = 0.0;
    double  userDt = solver.dt;
    double  outDt = solver.outDt;
    double  ftime = solver.ftime;
    double  dt;  /* delta time */

    /* save initial condition */
    swe_save_var(tstep++, time);
    phys->limit(phys, 2.0);
    swe_ppreserve(phys);
    double mpitime0 = MPI_Wtime();

    /* time step loop  */
    while (time<ftime){
        /* save result */
        if(time > outDt*tstep) {swe_save_var(tstep++, time);}
        /* calculate time interval */
        dt = swe_time_interval(phys);
        if(dt<userDt) dt = userDt;
        /* adjust final step to end exactly at FinalTime */
        if (time+dt > ftime) { dt = ftime-time; }
        if(!procid){ printf("Process:%f, dt:%f\n", time/ftime, dt); }

        for (INTRK=1; INTRK<=5; ++INTRK) {
            /* compute rhs of equations */
            const dg_real fdt = (dg_real)dt;
            const dg_real fa = (dg_real)rk4a[INTRK-1];
            const dg_real fb = (dg_real)rk4b[INTRK-1];

            phys->obc_update(phys, time+rk4c[INTRK-1]*dt);
            swe_rhs(phys, fa, fb, fdt);
            //swe_h2eta(solver);
            phys->limit(phys, 1.0);
            //swe_eta2h(solver);
            swe_ppreserve(phys);
        }
        /* increment current time */
        time += dt;
    }
    double mpitime1 = MPI_Wtime();
    double elapsetime = mpitime1 - mpitime0;

    // last
    swe_save_var(tstep++, time);

    if(!procid) {printf("proc: %d, time taken: %lg\n", procid, elapsetime);}
    return;
}

static void swe_rk_parameter(double *rk4a, double *rk4b, double *rk4c){
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
 * @brief Predict the wave speed and give local delta time
 *
 * @details
 * The wave speed is given as \f[ c = \sqrt{gh}+\left| \mathbf{u} \right| \f],
 * and the delta time is derived \f[ dt = dx/dt \f], while \f$ dx \f$ is the
 * characteristic length.
 *
 * @param[in] phys pointer to dg_phys structure;
 *
 * @return dt calculation time interval.
 */
static double swe_time_interval(dg_phys *phys){

    extern SWE_Solver solver;

    double dt = 1e4;
    const double gra  = solver.gra;
    const double hcrit = solver.hcrit;

    const int Np = dg_cell_Np( dg_phys_cell(phys) );
    const int K = dg_grid_K( dg_phys_grid(phys) );
    const int Nfield = dg_phys_Nfield(phys);
    dg_real *f_Q = dg_phys_f_Q(phys);
    register int k,n;
    dg_real h,u,v;

    for(k=0;k<K;k++){
        double len = dg_region_len(dg_phys_region(phys))[k];
        double c = 1e-10;
        int isdry = 0;
        for(n=0;n<Np;n++){
            int ind = (k*Np+n)*Nfield;
            h = f_Q[ind++];
            u = f_Q[ind++]; u /= h;
            v = f_Q[ind  ]; v /= h;

            if(h<hcrit){ isdry = 1; }
        c = max(c, sqrt(gra*h)+sqrt(u*u+v*v));
    }
    if(isdry){ continue; } // jump this cell
        dt = min(dt, len/c);
#if 0
        int procid;
        MPI_Comm_rank(MPI_COMM_WORLD, &procid);
        if(!procid) printf("k=%d, hmean=%f, c=%f, dt=%f\n",k,hmean,c,dt);
#endif
    }
    /* gather all the delta time in all process */
    double gdt;
    dt = dt*solver.cfl;
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
 * @param[in,out] phys pointer to dg_phys structure;
 */
static void swe_ppreserve(dg_phys *phys){

    const int Nfield = dg_phys_Nfield(phys);
    const int Np = dg_cell_Np( dg_phys_cell(phys) );
    const int K = dg_grid_K( dg_phys_grid(phys) );

    extern SWE_Solver solver;
    const double hcrit = solver.hcrit;

    register int i,k,ind;

    phys->cell_mean(phys);
    dg_real *c_Q = dg_phys_c_Q(phys);
    dg_real *f_Q = dg_phys_f_Q(phys);
    for(k=0;k<K;k++){
        ind = k*Nfield; /* index of k-th cell */
        dg_real hmean  = c_Q[ind];
        dg_real qxmean = c_Q[ind+1];
        dg_real qymean = c_Q[ind+2];
        /* correct negative water depth */
        if(hmean<0.0){
            for(i=0;i<Np;i++) {
                ind = (k*Np + i)*Nfield;
                f_Q[ind  ] -= hmean;
                f_Q[ind+1]  = 0.0;
                f_Q[ind+2]  = 0.0;
            }
            hmean = 0.0;
            qxmean = 0.0;
            qymean = 0.0;
        }

        ind = k*Np*Nfield; /* index of first node */
        dg_real hmin = f_Q[ind];
        /* compute minimum water depth */
        for(i=1;i<Np;i++){
            ind += Nfield;
            hmin = min(hmin, f_Q[ind]);
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
            f_Q[ind  ] = (f_Q[ind  ] - hmean )*theta + hmean;
            f_Q[ind+1] = (f_Q[ind+1] - qxmean)*theta + qxmean;
            f_Q[ind+2] = (f_Q[ind+2] - qymean)*theta + qymean;

            if(f_Q[ind]<hcrit){
                f_Q[ind+1] = 0.0;
                f_Q[ind+2] = 0.0;
            }
        }
    }

    return;
}
