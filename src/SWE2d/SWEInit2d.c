#include "SWEDriver2d.h"

void ParabolicBowlInit(SWESolver *solver, PhysDomain2d *phys, MultiReg2d *mesh);
void DamBreakDryInit  (SWESolver *solver, PhysDomain2d *phys, MultiReg2d *mesh);
void DamBreakWetInit  (SWESolver *solver, PhysDomain2d *phys, MultiReg2d *mesh);

/**
 * @brief
 * Initialize the variables and parameters for test cases.
 *
 * @details
 * Initialize the variables and set the final time for the simulation.
 *
 *
 * @param[in] casename
 * @param[in] solver
 * @param[in] mesh
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * phys | PhysDomain2d |
 *
 */
PhysDomain2d* SWEInit2d(char *casename, SWESolver *solver, MultiReg2d *mesh){
    PhysDomain2d *phys = GetPhysDomain2d(mesh, 3); /* 3 variables */

    /* initialize variables for various tests */
    if      ( !(memcmp(casename, "ParabolicBowl", 13)) ){
        ParabolicBowlInit(solver, phys, mesh);
    }else if( !(memcmp(casename, "DamBreakDry"  , 11)) ){
        DamBreakDryInit(solver, phys, mesh);
    }else if( !(memcmp(casename, "DamBreakWet"  , 11)) ){
        DamBreakWetInit(solver, phys, mesh);
    }else{
        printf("Wrong name of test case: %s\n", casename);
        MPI_Finalize(); exit(1);
    }

    return phys;
}

void ParabolicBowlInit(SWESolver *solver, PhysDomain2d *phys, MultiReg2d *mesh){
    double alpha = 1.6e-7;
    double w     = sqrt(8.0*solver->gra*alpha);
    double T     = 2*M_PI/w;
    double X     = 1.0;
    double Y     = -0.41884;

    StdRegions2d *shape = mesh->stdcell;

    /* finial time */
    solver->FinalTime = T;

    /* initialization */
    int k, i, ind, Nfields = 3;
    double r2, x2, y2;
    for (k=0; k<mesh->K; k++){
        for (i=0; i<shape->Np; i++){
            x2 = mesh->x[k][i]*mesh->x[k][i];
            y2 = mesh->y[k][i]*mesh->y[k][i];
            r2 = (X+Y)/(alpha*(X*X - Y*Y));
            ind = (k*shape->Np + i)*Nfields;
            // printf("k = %d, i = %d, ind = %d\n", k, i, ind);
            if( (x2 + y2) > r2 ){
                phys->f_Q[ind++] = 0;
            }else{
                phys->f_Q[ind++] = (real)( 1.0/(X+Y) + alpha*(Y*Y - X*X)*(x2 + y2)/(X+Y)/(X+Y) );
            }

            phys->f_Q[ind++] = 0; /* variable for qx */
            phys->f_Q[ind++] = 0; /* variable for qy */
        }
    }
}

void DamBreakDryInit(SWESolver *solver, PhysDomain2d *phys, MultiReg2d *mesh){
    StdRegions2d *shape = mesh->stdcell;
    double damPosition  = 500.0;

    /* finial time */
    solver->FinalTime = 20.0;

    /* initialization */
    int i,k,j,ind;
    int Nfields = 3;
    double xc = 0.0;
    for (k=0; k<mesh->K; k++){
        /* element centre */
        xc = 0;
        for (j=0;j<shape->Np;j++){
            xc += mesh->x[k][j];
        }
        xc /= (double)shape->Np;

        for (i=0; i<shape->Np; i++){
            ind = (k*shape->Np + i)*Nfields;

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


void DamBreakWetInit(SWESolver *solver, PhysDomain2d *phys, MultiReg2d *mesh){
    StdRegions2d *shape = mesh->stdcell;
    double damPosition  = 500.0;

    /* finial time */
    solver->FinalTime = 20.0;

    /* initialization */
    int i,k,j,ind;
    int Nfields = 3;
    double xc = 0;
    for (k=0; k<mesh->K; k++){

        /* element centre */
        xc = 0;
        for (j=0;j<shape->Np;j++){
            xc += mesh->x[k][j];
        }
        xc /= (double)shape->Np;

        for (i=0; i<shape->Np; i++){

            ind = (k*shape->Np + i)*Nfields;
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