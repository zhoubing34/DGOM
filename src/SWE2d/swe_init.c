#include "swe_init.h"
#include "swe_read_command.h"
#include "swe_mesh.h"
#include "swe_input.h"
#include "PhysField/dg_phys_init.h"

#define DEBUG 0
static void swe_initial_condition(swe_solver *solver);
static void swe_dambreakdry_init(swe_solver *solver);
static void swe_dambreakwet_init(swe_solver *solver);
static void swe_parabolicbowl_init(swe_solver *solver);
static void swe_userset_init(swe_solver *solver);
static void swe_set_topography(swe_solver *solver);

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
            swe_set_topography(solver);
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
        case swe_userset:
            swe_userset_init(solver); break;
    }
}

static void swe_set_topography(swe_solver *solver){
    char bot_filename[MAX_NAME_LENGTH];
    strcpy(bot_filename, solver->casename);
    strcat(bot_filename, ".bot");
    switch (solver->caseid){
        case swe_dambreakwet:
            solver->bot = swe_flat_topography(solver); break;
        case swe_dambreakdry:
            solver->bot = swe_flat_topography(solver); break;
        case swe_parabolicbowl:
            solver->bot = swe_parabolic_topography(solver); break;
        case swe_userset:
            solver->bot = swe_read_topography(solver, bot_filename); break;
    }
}

static void swe_userset_init(swe_solver *solver){
    physField *phys = solver->phys;
    pf_init_file2d(phys, solver->casename);
}

static void swe_dambreakdry_init(swe_solver *solver){

    physField *phys = solver->phys;
    const int K = phys->grid->K;
    const int Np = phys->cell->Np;

    const double x0 = 500; // dam position
    register int k,n,sk=0;

    dg_real *f_Q = phys->f_Q;
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

    dg_real *f_Q = phys->f_Q;
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

    dg_real *f_Q = phys->f_Q;
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