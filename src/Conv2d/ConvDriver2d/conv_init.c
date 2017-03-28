//
// Created by li12242 on 16/12/27.
//
#include "conv2d.h"

#define DEBUG 0

static dg_grid* user_grid_init();
static dg_phys* user_phys_init(dg_grid *grid);
static void conv_time_init(dg_phys *phys);

typedef struct Conv_Init_Creator{
    dg_grid *(*grid_init)();
    dg_phys *(*phys_init)(dg_grid *grid);
    void (*time_init)(dg_phys *phys);
}Conv_Init_Creator;

static const Conv_Init_Creator user_set_creator = {
        user_grid_init,
        user_phys_init,
        conv_time_init,
};

void conv_init(){
    extern Conv_Solver solver;
    const Conv_Init_Creator *solver_creator = &user_set_creator;

    dg_grid *grid = solver_creator->grid_init();
    dg_phys *phys = solver_creator->phys_init(grid);
    solver_creator->time_init(phys);
    solver.phys = phys;
    return;
}

static dg_grid* user_grid_init(){
    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    extern Conv_Solver solver;
    // read input file
    arg_section **sec_p = conv_read_inputfile(solver.filename);
    /// 0. case name
    char casename[MAX_NAME_LENGTH];
    arg_section *sec = sec_p[0];
    strcpy(casename, sec->arg_vec_p[0]);
    if(!procid){ printf(HEADLINE " casename: %s\n", casename); }
    /// 1. cell info
    sec = sec_p[1];
    dg_cell_type cell_type;
    int N;
    sscanf(sec->arg_vec_p[0], "%d\n", &(cell_type));
    sscanf(sec->arg_vec_p[1], "%d\n", &(N));
    dg_cell *cell = dg_cell_creat(N, cell_type);
    if(!procid){
        switch (cell_type){
            case TRIANGLE:
                printf(HEADLINE " cell type: triangle\n"); break;
            case QUADRIL:
                printf(HEADLINE " cell type: quadrilateral\n"); break;
            default:
                fprintf(stderr, "%s (%d)\nUnknown cell type %d.\n",
                        __FUNCTION__, __LINE__, cell_type);
                MPI_Abort(MPI_COMM_WORLD, -1);
        }
        printf(HEADLINE " polynomial degree: %d\n", N);
    }
    conv_arg_section_free(sec_p);

    // create grid
    dg_grid *grid = dg_grid_read_file2d(cell, casename);
    dg_grid_add_BS_file2d(grid, casename); // read boundary condition

    return grid;
}

static dg_phys* user_phys_init(dg_grid *grid){
    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    dg_region *region = dg_region_create(grid);
    dg_mesh *mesh = dg_mesh_create(region);
    dg_edge *edge = dg_edge_create(mesh);
    dg_phys *phys = dg_phys_create(3, edge);
    /// 2. obc file
    extern Conv_Solver solver;
    arg_section **sec_p = conv_read_inputfile(solver.filename);
    char obcfile[MAX_NAME_LENGTH];
    arg_section *sec = sec_p[2];
    strcpy(obcfile, sec->arg_vec_p[0]);
    if(!procid) printf(HEADLINE " open boundary file: %s\n", obcfile);
    if(strlen(obcfile)) { phys->obc_add(phys, obcfile); }

    /// 4. initial condition
    char filename[MAX_NAME_LENGTH];
    sec = sec_p[4];
    strcpy(filename, sec->arg_vec_p[0]);
    if(!procid) printf(HEADLINE " initial condition file: %s\n", filename);
    phys->init_file(phys, filename);

    conv_arg_section_free(sec_p);
    return phys;
}

static void conv_time_init(dg_phys *phys){
    double dt = 1e6; // initial time step
    const int N = dg_cell_N(dg_phys_cell(phys));
    const int Np = dg_cell_Np(dg_phys_cell(phys));
    const int K = dg_grid_K(dg_phys_grid(phys));
    double *len = dg_phys_region(phys)->len;

    dg_real *f_Q = dg_phys_f_Q(phys);
    int k,n,sk = 0;
    for(k=0;k<K;++k){
        double r = len[k]/(N+1);
        for(n=0;n<Np;n++){
            sk++; // jump c field
            const dg_real u = f_Q[sk++];
            const dg_real v = f_Q[sk++];
            double spe = sqrt(u*u+v*v);
            dt = min(dt, r/spe);
        }
    }
    /// 3. time info, read CFL number and final time
    extern Conv_Solver solver;
    arg_section **sec_p = conv_read_inputfile(solver.filename);

    double cfl,dt_user,ftime;
    arg_section *sec = sec_p[3];
    sscanf(sec->arg_vec_p[0], "%lf\n", &(cfl));
    sscanf(sec->arg_vec_p[1], "%lf\n", &(dt_user));
    sscanf(sec->arg_vec_p[2], "%lf\n", &(ftime));
    conv_arg_section_free(sec_p);

    if(dt_user < dt) {dt = dt_user;}
    double gdt;
    MPI_Allreduce(&dt, &gdt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt = cfl*gdt;

    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    if(!procid) printf(HEADLINE " dt: %f\n", dt);
    if(!procid) printf(HEADLINE " final time: %f\n", ftime);
    // assignment
    extern Conv_Solver solver;
    solver.finaltime = ftime;
    solver.dt = dt;
    return;
}
