#include "SWEDriver2d.h"

void ParabolicBowlInit2d(SWE_Solver2d *solver, PhysDomain2d *phys, MultiReg2d *mesh, double hmin);
void DamBreakDryInit2d(SWE_Solver2d *solver, PhysDomain2d *phys, MultiReg2d *mesh, double hmin);
void DamBreakWetInit2d(SWE_Solver2d *solver, PhysDomain2d *phys, MultiReg2d *mesh, double hmin);

/**
 * @brief
 * Initialize the variables and parameters for different test cases.
 *
 * @details
 * The parameters including:
 *  1. the critical depth for determing wet or dry;
 *  2. acceleration of gravity;
 *  3. minimum time step;
 *  4. simulation time.
 *
 * @param[in] argv input argument
 * @param[in] solver  SWE solver object
 * @param[in] mesh    mesh object.
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * phys | PhysDomain2d | pointer to physical domain object
 *
 */
PhysDomain2d* SWE_Init2d(char **argv, SWE_Solver2d *solver, MultiReg2d *mesh){
    PhysDomain2d *phys = GenPhysDomain2d(mesh, 3); /* 3 variables */

    double hmin = 1.0e-3; /* minimum depth */
    /* initialize variables and parameters for various tests */
    if      ( !(memcmp(argv[1], "ParabolicBowl", 13)) ){
        ParabolicBowlInit2d(solver, phys, mesh, hmin);
    }else if( !(memcmp(argv[1], "DamBreakDry"  , 11)) ){
        DamBreakDryInit2d(solver, phys, mesh, hmin);
    }else if( !(memcmp(argv[1], "DamBreakWet"  , 11)) ){
        DamBreakWetInit2d(solver, phys, mesh, hmin);
    }else{
        printf("Wrong name of test case: %s\n", argv[1]);
        MPI_Finalize(); exit(1);
    }
    return phys;
}

void ParabolicBowlInit2d(SWE_Solver2d *solver, PhysDomain2d *phys, MultiReg2d *mesh, double hmin){
    double alpha = 1.6e-7;
    double w     = sqrt(8.0*solver->gra*alpha);
    double T     = 2*M_PI/w;
    double X     = 1.0;
    double Y     = -0.41884;

    StdRegions2d *shape = mesh->stdcell;

    /* finial time */
    solver->FinalTime = T;
    /* minimum depth */
    solver->hcrit = hmin;
    /* minimum time step */
    solver->dtmin = 1.0e-4;

    /* initialization of the variables */
    int k, i, ind;
    double x2, y2;
    double r2 = (X+Y)/(alpha*(X*X - Y*Y));
    for (ind=0,k=0; k<mesh->K; k++){
        for (i=0; i<shape->Np; i++){
            x2 = mesh->x[k][i]*mesh->x[k][i];
            y2 = mesh->y[k][i]*mesh->y[k][i];
            // printf("k = %d, i = %d, ind = %d\n", k, i, ind);
            if( (x2 + y2) > r2 ){
                phys->f_Q[ind++] = 1e-6;
            }else{
                phys->f_Q[ind++] = (real)( 1.0/(X+Y) + alpha*(Y*Y - X*X)*(x2 + y2)/(X+Y)/(X+Y) );
            }
            phys->f_Q[ind++] = 0; /* variable for qx */
            phys->f_Q[ind++] = 0; /* variable for qy */
        }
    }
}

void DamBreakDryInit2d(SWE_Solver2d *solver, PhysDomain2d *phys, MultiReg2d *mesh, double hmin){
    StdRegions2d *shape = mesh->stdcell;
    double damPosition  = 500.0;

    /* finial time */
    solver->FinalTime = 20.0;
    /* minimum depth */
    solver->hcrit = hmin;
    /* minimum time step */
    solver->dtmin = 1.0e-6;

    /* initialization */
    int i,k,ind;
    double xc;
    for (ind=0,k=0; k<mesh->K; k++){
        /* element centre */
        xc = 0;
        for (i=0;i<shape->Np;i++){
            xc += mesh->x[k][i];
        }
        xc /= (double)shape->Np;

        for (i=0; i<shape->Np; i++){
            if(xc < damPosition){
                phys->f_Q[ind++] = 10.0;
            }else{
                phys->f_Q[ind++] = 0.0;
            }
            phys->f_Q[ind++] = 0; /* variable for qx */
            phys->f_Q[ind++] = 0; /* variable for qy */
        }
    }
}


void DamBreakWetInit2d(SWE_Solver2d *solver, PhysDomain2d *phys, MultiReg2d *mesh, double hmin){
    StdRegions2d *shape = mesh->stdcell;
    double damPosition  = 500.0;

    /* finial time */
    solver->FinalTime = 20.0;
    /* minimum depth */
    solver->hcrit = hmin;
    /* minimum time step */
    solver->dtmin = 1.0e-4;

    /* initialization */
    int i,k,ind;
    double xc;
    for (ind=0,k=0; k<mesh->K; k++){
        /* element centre */
        xc = 0;
        for (i=0;i<shape->Np;i++){
            xc += mesh->x[k][i];
        }
        xc /= (double)shape->Np;

        for (i=0; i<shape->Np; i++){

            if(xc < damPosition){
                phys->f_Q[ind++] = 10.0;
            }else{
                phys->f_Q[ind++] = 2.0;
            }

            phys->f_Q[ind++] = 0; /* variable for qx */
            phys->f_Q[ind++] = 0; /* variable for qy */
        }
    }
}