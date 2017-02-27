//
// Created by li12242 on 17/1/6.
//

#include <PhysField/pf_phys.h>
#include <MultiRegions/mr_grid.h>
#include "conv_extsol.h"

extern conv_solver2d solver;

/**
 * @brief get exact solution in each element
 * @param [in] Np number of points in each
 * @param [in] x coordinate
 * @param [in] y coordinate
 * @param [in,out] ext exact solution
 */
static void advection_diffusion_ext(int Np, dg_real *x, dg_real *y, double *ext){

    const double t = solver.finaltime;
    const double u = solver.u, v = solver.v;
    const double x0 = -0.5, y0 = -0.5;

    int register n;
    if(solver.viscosity>DIFF_THRESHOLD){
        for(n=0;n<Np;n++){
            const double miu = 1.0/solver.viscosity;
            const double tiff = 1.0/(4.0*t+1);
            double xiff = x[n]-x0-u*t;
            double yiff = y[n]-y0-v*t;
            double cx = sqrt(tiff)*exp(-tiff*xiff*xiff*miu);
            double cy = sqrt(tiff)*exp(-tiff*yiff*yiff*miu);
            ext[n] = cx*cy;
        }
    }else{
        const double sigma = 125*1e3/(33*33);
        for(n=0;n<Np;n++){
            double xiff = x[n]-x0-u*t;
            double yiff = y[n]-y0-v*t;
            double cx = exp(-sigma*xiff*xiff);
            double cy = exp(-sigma*yiff*yiff);
            ext[n] = cx*cy;
        }
    }


    return;
}

/**
 * @brief exact solution in each element for rotation case
 * @param [in] Np number of points
 * @param [in] x coordinate
 * @param [in] y coordinate
 * @param [in,out] ext exact solution
 */
static void rotation_ext(int Np, dg_real *x, dg_real *y, double *ext){

    const double t = solver.finaltime;
    const double r = 0.6, phase = M_PI_2, T = 2.4;
    const double theta = phase + t/T*2*M_PI;
    const double x0 = cos(theta)*r;
    const double y0 = sin(theta)*r;

    const double sigma = 125*1e3/(33*33);

    int register n;
    for(n=0;n<Np;n++){
        const double xt = x[n];
        const double yt = y[n];
        double c = -sigma * ( ( xt - x0 )*( xt - x0 ) + ( yt - y0 )*( yt - y0 ) );
        ext[n] = exp(c);
    }
    return;
}


typedef void (*extsol_fun)(int Np, dg_real *x, dg_real *y, double *ext);
/**
 * @brief
 * @details
 */
void conv_normerr(physField *phys){

    const int K = phys->grid->K;
    const int Nfield = phys->Nfield;
    const int Np = phys->cell->Np;
    const int procid = phys->grid->procid;

    double Linf=0, L2=0, L1=0;

    extsol_fun extsolFun=NULL;
    switch (solver.caseid){
        case conv_rotational_convection:
            extsolFun = rotation_ext; break;
        case conv_advection_diffusion:
            extsolFun = advection_diffusion_ext; break;
        default:
            if(!procid) { printf("%s: %d\nUnknown exact solution, exit\n", __FUNCTION__, __LINE__);}
            return;
    }

    int register k,n,sk=0;
    dg_real *f_Q = phys->f_Q;
    double var[Np], ext_Q[Np];
    double err1[Np], err2[Np];
    double Aloc=0; /* total area */

    for(k=0;k<K;k++){
        extsolFun(Np, phys->region->x[k], phys->region->y[k], ext_Q);
        for(n=0;n<Np;n++){
            var[n] = (double) f_Q[sk]; sk+=Nfield;

            err1[n] = fabs(var[n] - ext_Q[n]);
            err2[n] = err1[n]*err1[n];
            Linf = max(Linf, err1[n]);
        }

        Aloc += phys->region->size[k];

        L1 += mr_reg_integral(phys->region, k, err1);
        L2 += mr_reg_integral(phys->region, k, err2);
    }
    double gL1, gLinf, gL2, Atol;

    MPI_Allreduce(&Linf, &gLinf, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&L2, &gL2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&L1, &gL1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&Aloc, &Atol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    gL2 = sqrt(gL2/Atol);
    gL1 /= Atol;

    if(!procid)
        printf("proc: %d,\t L1: %lg,\t L2: %lg,\t Linf: %lg\n",
               procid, gL1, gL2, gLinf);

    return;
}

