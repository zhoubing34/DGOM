//
// Created by li12242 on 17/1/18.
//

#include "swe_extsol.h"

typedef void (*extsol_fun)(dg_real x, dg_real y, double t, dg_real *ext);

static double gra;

void swe_parabolicbowl_ext(dg_real x, dg_real y, double t, dg_real *ext){

    //extern double gra;
    const double alpha = 1.6*1e-7;
    const double w = sqrt(8*gra*alpha);
    const double X = -1;
    const double Y = -0.41884;

    const double r2 = x*x + y*y;
    const double r2ext = (X+Y*cos(w*t))/(alpha*(X*X - Y*Y));
    if(r2 > r2ext){
        ext[0] = 0;
        ext[1] = 0;
        ext[2] = 0;
    }else{
        double temp = (X+Y*cos(w*t))*(X+Y*cos(w*t));
        ext[0] = 1/(X+Y*cos(w*t)) + alpha*(Y*Y - X*X)*r2/temp;
        temp = -(Y*w*sin(w*t))/(X+Y*cos(w*t))*x/2.0;
        ext[1] = temp*ext[0];
        temp = -(Y*w*sin(w*t))/(X+Y*cos(w*t))*y/2.0;
        ext[2] = temp*ext[0];
    }
}

void swe_normerr(SWE_Solver *solver){

    const physField *phys = solver->phys;
    const int K = phys->grid->K;
    const int Nfield = phys->Nfield;
    const int Np = phys->cell->Np;

    /* set global variable */
    //extern double gra;
    gra = solver->gra;

    extsol_fun extsolFun=NULL;

    switch (solver->caseid){
        case swe_dambreakdry:
            break;
        case swe_dambreakwet:
            break;
        case swe_parabolicbowl:
            extsolFun = swe_parabolicbowl_ext; break;
        default:
            return;
    }

    int register k,n,fld,sk=0;
    const double ftime = solver->ftime;
    dg_real *f_Q = phys->f_Q;
    double ext_Q[Nfield], var, Aloc=0;
    double **err1 = matrix_double_create(Nfield, Np);
    double **err2 = matrix_double_create(Nfield, Np);

    double *Linf = vector_double_create(Nfield);
    double *L1 = vector_double_create(Nfield);
    double *L2 = vector_double_create(Nfield);

    for(k=0;k<K;k++){
        double *x = phys->region->x[k];
        double *y = phys->region->y[k];
        for(n=0;n<Np;n++){
            extsolFun(x[n], y[n], ftime, ext_Q);
            for(fld=0;fld<Nfield;fld++){
                var = (double) f_Q[sk++];
                double temp = fabs(var - ext_Q[fld]);
                err1[fld][n] = temp;
                err2[fld][n] = temp*temp;
                Linf[fld] = max(Linf[fld], temp);
            }
        }

        Aloc += phys->region->size[k];

        for(fld=0;fld<Nfield;fld++){
            L1[fld] += mr_reg_integral(phys->region, k, err1[fld]);
            L2[fld] += mr_reg_integral(phys->region, k, err2[fld]);
        }

    }
    const int procid = phys->mesh->procid;
    double Atol;
    MPI_Allreduce(&Aloc, &Atol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double gL1[Nfield], gLinf[Nfield], gL2[Nfield];

    for(fld=0;fld<Nfield;fld++) {
        MPI_Allreduce(Linf+fld, gLinf+fld, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(L2+fld, gL2+fld, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(L1+fld, gL1+fld, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        gL2[fld] = sqrt(gL2[fld] / Atol);
        gL1[fld] /= Atol;
    }
    if (!procid){
        printf("h: L1\t L2\t Linf\n%lg  %lg  %lg\n", gL1[0], gL2[0], gLinf[0]);
        printf("qx: L1\t L2\t Linf\n%lg  %lg  %lg\n", gL1[1], gL2[1], gLinf[1]);
        printf("qy: L1\t L2\t Linf\n%lg  %lg  %lg\n", gL1[2], gL2[2], gLinf[2]);
    }

    vector_double_free(L1);
    vector_double_free(L2);
    vector_double_free(Linf);

    matrix_double_free(err1);
    matrix_double_free(err2);
    return;
}