#include "swe_init.h"
#include "swe_read_command.h"
#include "swe_input.h"
#include "swe_driver2d.h"

#define DEBUG 0
static void swe_initial_condition(swe_solver *solver);
static void swe_dambreakdry_init(swe_solver *solver);
static void swe_dambreakwet_init(swe_solver *solver);
static void swe_parabolicbowl_init(swe_solver *solver);

swe_solver* swe_init(int argc, char **argv){
    swe_command_type type = swe_read_command(argc, argv);

    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    swe_solver *solver = NULL;
    switch (type) {
        case swe_command_help:
            if(!procid) {swe_print_help();} exit(0); // print help info and exit
        case swe_command_create_input:
            if(!procid) {swe_create_input();} exit(0); // create input file and exit
        case swe_command_run:
            solver = swe_create_solver();
            swe_initial_condition(solver);
            break;
        default:
            fprintf(stderr, "Unknown run type %d\n", type);
            MPI_Abort(MPI_COMM_WORLD, -1);
    }
    return solver;
}

/**
 * @brief create initial condition
 */
static void swe_initial_condition(swe_solver *solver){
    switch (solver->caseid){
        case swe_dambreakdry:
            swe_dambreakdry_init(solver); break;
        case swe_dambreakwet:
            swe_dambreakwet_init(solver); break;
        case swe_parabolicbowl:
            swe_parabolicbowl_init(solver); break;
    }

}

static void swe_dambreakdry_init(swe_solver *solver){

    physField *phys = solver->phys;
    const int K = phys->grid->K;
    const int Np = phys->cell->Np;

    const double x0 = 500; // dam position
    register int k,n,sk=0;

    real *f_Q = phys->f_Q;
    for(k=0;k<K;k++){
        double *xt = phys->region->x[k];
        double area = phys->region->size[k];
        double xc = mr_reg_integral(phys->region, k, xt)/area;

        for(n=0;n<Np;n++){
            // h field
            if(xc<x0){
                f_Q[sk++] = 10.0;
            }else{
                f_Q[sk++] = 0.0;
            }
            f_Q[sk++] = 0.0; // q_x field
            f_Q[sk++] = 0.0; // q_y field
        }
    }
    return;
}

static void swe_dambreakwet_init(swe_solver *solver){

    physField *phys = solver->phys;
    const int K = phys->grid->K;
    const int Np = phys->cell->Np;

    const double x0 = 500; // dam position
    register int k,n,sk=0;

    real *f_Q = phys->f_Q;
    for(k=0;k<K;k++){
        double *xt = phys->region->x[k];
        double area = phys->region->size[k];
        double xc = mr_reg_integral(phys->region, k, xt)/area;

        for(n=0;n<Np;n++){
            // h field
            if(xc<x0){
                f_Q[sk++] = 10.0;
            }else{
                f_Q[sk++] = 2.0;
            }
            f_Q[sk++] = 0.0; // q_x field
            f_Q[sk++] = 0.0; // q_y field
        }
    }
    return;
}

static void swe_parabolicbowl_init(swe_solver *solver){

    physField *phys = solver->phys;
    const int K = phys->grid->K;
    const int Np = phys->cell->Np;

    const double alpha = 1.6e-7;
    const double X = 1;
    const double Y = -0.41884;

    real *f_Q = phys->f_Q;
    register int k,n,sk=0;
    for(k=0;k<K;k++){
        for(n=0;n<Np;n++){
            const double xt = phys->region->x[k][n];
            const double yt = phys->region->y[k][n];
            const double r2 = (xt*xt + yt*yt);
            const double r2ext = (X+Y)/(alpha*(X*X - Y*Y));
            if(r2<r2ext){
                f_Q[sk++] = 1.0/(X+Y) + alpha*(Y*Y - X*X)*r2/(X+Y)/(X+Y);
                f_Q[sk++] = 0.0;
                f_Q[sk++] = 0.0;
            }else{
                f_Q[sk++] = 0.0;
                f_Q[sk++] = 0.0;
                f_Q[sk++] = 0.0;
            }
        }
    }
    return;
}

//void ParabolicBowlInit2d(SWE_Solver2d *solver, PhysDomain2d *phys, MultiReg2d *mesh, double hmin);
//void DamBreakDryInit2d(SWE_Solver2d *solver, PhysDomain2d *phys, MultiReg2d *mesh, double hmin);
//void DamBreakWetInit2d(SWE_Solver2d *solver, PhysDomain2d *phys, MultiReg2d *mesh, double hmin);
//
///**
// * @brief
// * Initialize the variables and parameters
// */
//PhysDomain2d* SWE_Init2d(char **argv, SWE_Solver2d *solver, MultiReg2d *mesh){
//    PhysDomain2d *phys = GenPhysDomain2d(mesh, 3); /* 3 variables */
//
//    double hmin = 1.0e-3; /* minimum depth */
//    /* initialize variables and parameters for various tests */
//    if      ( !(memcmp(argv[1], "ParabolicBowl", 13)) ){
//        ParabolicBowlInit2d(solver, phys, mesh, hmin);
//    }else if( !(memcmp(argv[1], "DamBreakDry"  , 11)) ){
//        DamBreakDryInit2d(solver, phys, mesh, hmin);
//    }else if( !(memcmp(argv[1], "DamBreakWet"  , 11)) ){
//        DamBreakWetInit2d(solver, phys, mesh, hmin);
//    }else{
//        printf("Wrong name of test case: %s\n", argv[1]);
//        MPI_Finalize(); exit(1);
//    }
//    return phys;
//}

//void ParabolicBowlInit2d(SWE_Solver2d *solver, PhysDomain2d *phys, MultiReg2d *mesh, double hmin){
//    double alpha = 1.6e-7;
//    double w     = sqrt(8.0*solver->gra*alpha);
//    double T     = 2*M_PI/w;
//    double X     = 1.0;
//    double Y     = -0.41884;
//
//    StdRegions2d *shape = mesh->stdcell;
//
//    /* finial time */
//    solver->FinalTime = T;
//    /* minimum depth */
//    solver->hcrit = hmin;
//    /* minimum time step */
//    solver->dtmin = 1.0e-4;
//
//    /* initialization of the variables */
//    int k, i, ind;
//    double x2, y2;
//    double r2 = (X+Y)/(alpha*(X*X - Y*Y));
//    for (ind=0,k=0; k<mesh->K; k++){
//        for (i=0; i<shape->Np; i++){
//            x2 = mesh->x[k][i]*mesh->x[k][i];
//            y2 = mesh->y[k][i]*mesh->y[k][i];
//            // printf("k = %d, i = %d, ind = %d\n", k, i, ind);
//            if( (x2 + y2) > r2 ){
//                phys->f_Q[ind++] = 1e-6;
//            }else{
//                phys->f_Q[ind++] = (real)( 1.0/(X+Y) + alpha*(Y*Y - X*X)*(x2 + y2)/(X+Y)/(X+Y) );
//            }
//            phys->f_Q[ind++] = 0; /* variable for qx */
//            phys->f_Q[ind++] = 0; /* variable for qy */
//        }
//    }
//}
//
//void DamBreakDryInit2d(SWE_Solver2d *solver, PhysDomain2d *phys, MultiReg2d *mesh, double hmin){
//    StdRegions2d *shape = mesh->stdcell;
//    double damPosition  = 500.0;
//
//    /* finial time */
//    solver->FinalTime = 20.0;
//    /* minimum depth */
//    solver->hcrit = hmin;
//    /* minimum time step */
//    solver->dtmin = 1.0e-6;
//
//    /* initialization */
//    int i,k,ind;
//    double xc;
//    for (ind=0,k=0; k<mesh->K; k++){
//        /* element centre */
//        xc = 0;
//        for (i=0;i<shape->Np;i++){
//            xc += mesh->x[k][i];
//        }
//        xc /= (double)shape->Np;
//
//        for (i=0; i<shape->Np; i++){
//            if(xc < damPosition){
//                phys->f_Q[ind++] = 10.0;
//            }else{
//                phys->f_Q[ind++] = 0.0;
//            }
//            phys->f_Q[ind++] = 0; /* variable for qx */
//            phys->f_Q[ind++] = 0; /* variable for qy */
//        }
//    }
//}
//
//
//void DamBreakWetInit2d(SWE_Solver2d *solver, PhysDomain2d *phys, MultiReg2d *mesh, double hmin){
//    StdRegions2d *shape = mesh->stdcell;
//    double damPosition  = 500.0;
//
//    /* finial time */
//    solver->FinalTime = 20.0;
//    /* minimum depth */
//    solver->hcrit = hmin;
//    /* minimum time step */
//    solver->dtmin = 1.0e-4;
//
//    /* initialization */
//    int i,k,ind;
//    double xc;
//    for (ind=0,k=0; k<mesh->K; k++){
//        /* element centre */
//        xc = 0;
//        for (i=0;i<shape->Np;i++){
//            xc += mesh->x[k][i];
//        }
//        xc /= (double)shape->Np;
//
//        for (i=0; i<shape->Np; i++){
//
//            if(xc < damPosition){
//                phys->f_Q[ind++] = 10.0;
//            }else{
//                phys->f_Q[ind++] = 2.0;
//            }
//
//            phys->f_Q[ind++] = 0; /* variable for qx */
//            phys->f_Q[ind++] = 0; /* variable for qy */
//        }
//    }
//}