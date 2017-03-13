//
// Created by li12242 on 16/12/27.
//
#include "conv_driver.h"

#define DEBUG 0

static dg_grid* uniform_grid_init();
static dg_grid* user_grid_init();
static dg_phys* rotation_phys_init(dg_grid *grid);
static dg_phys* advection_phys_init(dg_grid *grid);
static dg_phys* user_phys_init(dg_grid *grid);
static void conv_time_init(dg_phys *phys);

typedef struct Conv_Init_Creator{
    dg_grid *(*grid_init)();
    dg_phys *(*phys_init)(dg_grid *grid);
    void (*time_init)(dg_phys *phys);
}Conv_Init_Creator;

static const Conv_Init_Creator rotation_creator = {
        uniform_grid_init,
        rotation_phys_init,
        conv_time_init,
};
static const Conv_Init_Creator advection_creator = {
        uniform_grid_init,
        advection_phys_init,
        conv_time_init,
};
static const Conv_Init_Creator user_set_creator = {
        user_grid_init,
        user_phys_init,
        conv_time_init,
};

void conv_init(Conv_Case_Type case_type){
    extern Conv_Solver solver;
    const Conv_Init_Creator *solver_creator;
    switch (case_type){
        case conv_rotational_convection:
            solver_creator = &rotation_creator; break;
        case conv_advection_diffusion:
            solver_creator = &advection_creator; break;
        case conv_userset:
            solver_creator = &user_set_creator; break;
        default:
            fprintf(stderr, "%s: %d\nUnknown case type %d\n",
                    __FUNCTION__, __LINE__, case_type);
            exit(-1);
    }

    dg_grid *grid = solver_creator->grid_init();
    dg_phys *phys = solver_creator->phys_init(grid);
    solver_creator->time_init(phys);
    return;
}

static dg_grid* uniform_grid_init(){
    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    // read input file
    arg_section **sec_p = conv_read_input();
    /// 0. case info
    /// 1. cell info
    arg_section *sec = sec_p[1];
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
    /// 2. grid info
    sec = sec_p[2];
    int Mx, My;
    sscanf(sec->arg_vec_p[0], "%d\n", &(Mx));
    sscanf(sec->arg_vec_p[1], "%d\n", &(My));
    dg_grid *grid = NULL;
    if(!procid){
        printf(HEADLINE " Ne on x: %d\n", Mx);
        printf(HEADLINE " Ne on y: %d\n", My);
    }
    switch (cell_type){
        case TRIANGLE:
            grid = dg_grid_uniform_tri(cell, Mx, My, -1, 1, -1, 1, 0); break;
        case QUADRIL:
            grid = dg_grid_uniform_quad(cell, Mx, My, -1, 1, -1, 1); break;
        default:
            fprintf(stderr, "%s (%d)\nUnknown cell type %d.\n",
                    __FUNCTION__, __LINE__, cell_type);
            MPI_Abort(MPI_COMM_WORLD, -1);
    }
    conv_arg_section_free(sec_p);
    return grid;
}

static dg_grid* user_grid_init(){
    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    // read input file
    arg_section **sec_p = conv_read_input();
    /// 0. case info
    /// 1. cell info
    arg_section *sec = sec_p[1];
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
    /// 2. grid info
    sec = sec_p[2];
    char *casename = (char *)calloc( strlen(sec->arg_vec_p[2])+1, sizeof(char));
    strcpy(casename, sec->arg_vec_p[2]);
    if(!procid){
        printf(HEADLINE " casename: %s\n", casename);
    }
    dg_grid *grid = dg_grid_read_file2d(cell, casename);
    dg_grid_read_BSfile2d(grid, casename); // read boundary condition
    conv_arg_section_free(sec_p);
    free(casename);
    return grid;
}

static dg_phys* rotation_phys_init(dg_grid *grid){
    dg_region *region = dg_region_create(grid);
    dg_mesh *mesh = dg_mesh_create(region);
    dg_edge *edge = dg_edge_create(mesh);
    dg_phys *phys = dg_phys_create(3, edge);
    const int K = dg_grid_K(phys->grid);
    const int Np = dg_cell_Np(phys->cell);
    const double sigma = 125*1e3/(33*33);
    const double w = 5*M_PI/6;
    const double xc = 0.0;
    const double yc = 0.6;

    double **x = phys->region->x;
    double **y = phys->region->y;

    int k,n,sk = 0;
    for(k=0;k<K;++k){
        for(n=0;n<Np;++n){
            const double xt = x[k][n];
            const double yt = y[k][n];
            double t = -sigma * (( xt - xc )*( xt - xc ) + ( yt - yc )*( yt - yc ));
            phys->f_Q[sk++] = (dg_real) exp(t); // c field
            phys->f_Q[sk++] = (dg_real)(-w * yt); // flow rate at x-coordinate
            phys->f_Q[sk++] = (dg_real)( w * xt); // flow rate at y-coordinate
        }
    }
    return phys;
}

static dg_phys* advection_phys_init(dg_grid *grid){
    dg_region *region = dg_region_create(grid);
    dg_mesh *mesh = dg_mesh_create(region);
    dg_edge *edge = dg_edge_create(mesh);
    dg_phys *phys = dg_phys_create(3, edge);

    const int K = dg_grid_K(phys->grid);
    const int Np = dg_cell_Np(phys->cell);
    const double xc = -0.5;
    const double yc = -0.5;
    const double sigma = 125*1e3/(33*33);

    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    double **x = phys->region->x;
    double **y = phys->region->y;
    // read flow field
    arg_section **sec_p = conv_read_input();
    /// 4. phys info
    arg_section *sec = sec_p[4];
    double ut,vt;
    sscanf(sec->arg_vec_p[0], "%lf\n", &ut);
    sscanf(sec->arg_vec_p[1], "%lf\n", &vt);
    if(!procid){
        printf(HEADLINE " u: %lf\n", ut);
        printf(HEADLINE " v: %lf\n", vt);
    }
    conv_arg_section_free(sec_p);

    int k,n,sk = 0;
    for(k=0;k<K;++k){
        for(n=0;n<Np;++n){
            const double xt = x[k][n];
            const double yt = y[k][n];
            double t = -( ( xt - xc )*( xt - xc ) + ( yt - yc )*( yt - yc ) )*sigma;
            phys->f_Q[sk++] = (dg_real) exp(t); // c field
            phys->f_Q[sk++] = (dg_real) ut; // flow rate at x-coordinate
            phys->f_Q[sk++] = (dg_real) vt; // flow rate at y-coordinate
        }
    }
    return phys;
}

static dg_phys* user_phys_init(dg_grid *grid){
    dg_region *region = dg_region_create(grid);
    dg_mesh *mesh = dg_mesh_create(region);
    dg_edge *edge = dg_edge_create(mesh);
    dg_phys *phys = dg_phys_create(3, edge);

    return phys;
}

static void conv_time_init(dg_phys *phys){
    double dt = 1e6; // initial time step
    const int N = dg_cell_N(phys->cell);
    const int Np = dg_cell_Np(phys->cell);
    const int K = dg_grid_K(phys->grid);
    double *len = phys->region->len;

    int k,n,sk = 0;
    for(k=0;k<K;++k){
        double r = len[k]/(N+1);
        for(n=0;n<Np;n++){
            sk++; // jump c field
            const dg_real u = phys->f_Q[sk++];
            const dg_real v = phys->f_Q[sk++];
            double spe = sqrt(u*u+v*v);
            dt = min(dt, r/spe);
        }
    }
    // read CFL number and final time
    arg_section **sec_p = conv_read_input();
    /// 3. time info
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

    // assignment
    extern Conv_Solver solver;
    solver.phys = phys;
    solver.finaltime = ftime;
    solver.dt = dt;
    return;
}