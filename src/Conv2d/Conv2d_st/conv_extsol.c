//
// Created by li12242 on 17/1/6.
//

#include "conv2d_st.h"

extern Conv_Solver solver;

/**
 * @brief get exact solution in each element
 * @param [in] Np number of points in each
 * @param [in] x coordinate
 * @param [in] y coordinate
 * @param [in,out] ext exact solution
 */
static void advdiff_ext(int Np, dg_real *x, dg_real *y, double *ext){

    const double t = solver.finaltime;
    const double x0 = -0.5, y0 = -0.5;
    int register n;
    dg_phys *phys = solver.phys;
    dg_real *f_Q = dg_phys_f_Q(phys);
    dg_phys_LDG *ldg = phys->ldg;
    const double u = f_Q[1];
    const double v = f_Q[2];

    //const double sigma = 125*1e3/(33*33);
    const double miu = 1.0/(dg_phys_ldg_sqrt_miux(ldg)[0]);
    const double tiff = 1.0/(4.0*t+1);
    for(n=0;n<Np;n++){
        double xiff = x[n]-x0-u*t;
        double yiff = y[n]-y0-v*t;
        double cx = sqrt(tiff)*exp(-tiff*xiff*xiff*miu);
        double cy = sqrt(tiff)*exp(-tiff*yiff*yiff*miu);
        ext[n] = cx*cy;
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

static void adv_ext(int Np, dg_real *x, dg_real *y, double *ext){
    const double t = solver.finaltime;
    dg_phys *phys = solver.phys;
    dg_real *f_Q = dg_phys_f_Q(phys);
    const double x0 = -0.5;
    const double y0 = -0.5;
    const double u = f_Q[1];
    const double v = f_Q[2];

    const double sigma = 125*1e3/(33*33);

    int register n;
    for(n=0;n<Np;n++){
        const double xt = x[n]-x0-u*t;
        const double yt = y[n]-y0-v*t;
        double c = -sigma * ( xt*xt + yt*yt );
        ext[n] = exp(c);
    }
    return;
}

typedef void (*Ext_Fun)(int Np, dg_real *x, dg_real *y, double *ext);

/**
 * @brief
 * @details
 */
void conv_normerr(){

    dg_phys *phys = solver.phys;
    dg_region *region = dg_phys_region(phys);
    const int K = dg_grid_K(dg_phys_grid(phys));
    const int Nfield = dg_phys_Nfield(phys);
    const int Np = dg_cell_Np(dg_phys_cell(phys));
    const int procid = dg_grid_procid(dg_phys_grid(phys));

    double Linf=0, L2=0, L1=0;
    double **x = dg_region_x(region);
    double **y = dg_region_y(region);

    Ext_Fun ext_fun=NULL;
    // read case ID
    arg_section **sec_p = conv_read_inputfile(solver.filename);
    arg_section *sec = sec_p[0];
    Conv_St_Case case_type;
    sscanf(sec->arg_vec_p[0], "%d\n", &(case_type));
    conv_arg_section_free(sec_p);

    switch (case_type){
        case Conv_Rotation:
            ext_fun = rotation_ext; break;
        case Conv_AdvDiff:
            ext_fun = advdiff_ext; break;
        case Conv_Adv:
            ext_fun = adv_ext; break;
        default:
            if(!procid){ printf("%s (%d)\nUnknown exact solution.\n", __FUNCTION__, __LINE__); }
            return;
    }

    int register k,n,sk=0;
    dg_real *f_Q = dg_phys_f_Q(phys);
    dg_real var[Np], ext_Q[Np];
    dg_real err1[Np], err2[Np];
    dg_real Aloc=0; /* total area */

    for(k=0;k<K;k++){
        ext_fun(Np, x[k], y[k], ext_Q);
        for(n=0;n<Np;n++){
            var[n] = f_Q[sk]; sk+=Nfield;
            err1[n] = fabs(var[n] - ext_Q[n]);
            err2[n] = err1[n]*err1[n];
            Linf = max(Linf, err1[n]);
        }

        dg_real tmp;
        Aloc += region->size[k];
        region->vol_integral(region, 1, k, err1, &tmp);
        L1 += tmp;
        region->vol_integral(region, 1, k, err2, &tmp);
        L2 += tmp;
    }
    double gL1, gLinf, gL2, Atol;

    MPI_Allreduce(&Linf, &gLinf, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&L2, &gL2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&L1, &gL1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&Aloc, &Atol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    gL2 = sqrt(gL2/Atol);
    gL1 /= Atol;

    if(!procid)
        printf("proc: %d,\t L1: %lg,\t L2: %lg,\t Linf: %lg\n", procid, gL1, gL2, gLinf);

    return;
}