//
// Created by li12242 on 17/3/28.
//
#include "conv2d_st.h"
#include "../ConvLib/conv_lib2d.h"
#include "MultiRegions/Grid/dg_grid_BS.h"

#define DEBUG 0
#if DEBUG
#include "Utility/unit_test.h"
#endif

typedef struct Conv_Init_Creator{
    dg_grid *(*grid_init)();
    dg_phys *(*phys_init)(dg_grid *grid);
    void (*time_init)(dg_phys *phys);
}Conv_Init_Creator;

static void conv_time_init(dg_phys *phys);
static void set_uniform_obc(dg_grid *grid, int Mx, int My);
static dg_grid* uniform_grid_init();
static dg_phys* rotation_phys_init(dg_grid *grid);
static dg_phys* advection_phys_init(dg_grid *grid);
static dg_phys* pure_advection_phys_init(dg_grid *grid);
static void advection_time_init(dg_phys *phys);
static void advdiff_time_init(dg_phys *phys);

static const Conv_Init_Creator rotation_creator = {
        uniform_grid_init,
        rotation_phys_init,
        advection_time_init,
};
static const Conv_Init_Creator advdiff_creator = {
        uniform_grid_init,
        advection_phys_init,
        advdiff_time_init,
};
static const Conv_Init_Creator pure_advection = {
        uniform_grid_init,
        pure_advection_phys_init,
        advection_time_init,
};


void conv_init_st(){
    extern Conv_Solver solver;
    const Conv_Init_Creator *solver_creator;

    /* read case ID */
    arg_section **sec_p = conv_read_inputfile(solver.filename);
    arg_section *sec = sec_p[0];
    Conv_St_Case case_type;
    sscanf(sec->arg_vec_p[0], "%d\n", &(case_type));
    conv_arg_section_free(sec_p);

    switch (case_type){
        case Conv_Rotation:
            solver_creator = &rotation_creator; break;
        case Conv_AdvDiff:
            solver_creator = &advdiff_creator; break;
        case Conv_Adv:
            solver_creator = &pure_advection; break;
        default:
            fprintf(stderr, "%s: %d\nUnknown case type %d\n",
                    __FUNCTION__, __LINE__, case_type);
            exit(-1);
    }

    dg_grid *grid = solver_creator->grid_init();
    dg_phys *phys = solver_creator->phys_init(grid);
    solver_creator->time_init(phys);
    solver.phys = phys;
    return;
}

static dg_grid* uniform_grid_init(){
    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    extern Conv_Solver solver;
    // read input file
    arg_section **sec_p = conv_read_inputfile(solver.filename);
    /// 0. case info
    arg_section *sec = sec_p[0];
    Conv_St_Case case_type;
    sscanf(sec->arg_vec_p[0], "%d\n", &(case_type));
    if(!procid){
        switch (case_type){
            case Conv_Rotation:
                printf(HEADLINE " case name: rotational convection\n"); break;
            case Conv_AdvDiff:
                printf(HEADLINE " case name: advection diffusion\n"); break;
            case Conv_Adv:
                printf(HEADLINE " case name: pure advection\n"); break;
            default:
                fprintf(stderr, "%s: %d\nUnknown case type %d\n",
                        __FUNCTION__, __LINE__, case_type);
                exit(-1);
        }
    }
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
    set_uniform_obc(grid, Mx, My);
    return grid;
}

static void set_uniform_obc(dg_grid *grid, int Mx, int My){
    const int Nsurf = 2*(Mx+My);
    int **SFToV = matrix_int_create(Nsurf, 3);
    int n, sf=0;
    /* bottom - inflow */
    for(n=0;n<Mx;n++){
        SFToV[sf  ][0] = n;
        SFToV[sf  ][1] = n+1;
        SFToV[sf++][2] = FACE_OPENBS;
    }
    /* top - outflow */
    int sk = (Mx+1)*My;
    for(n=0;n<Mx;n++){
        SFToV[sf  ][0] = n + sk;
        SFToV[sf  ][1] = n+1 + sk;
        SFToV[sf++][2] = FACE_OPENBS;
    }
    /* left - inflow */
    sk = Mx+1;
    for(n=0;n<My;n++){
        SFToV[sf  ][0] = n*sk;
        SFToV[sf  ][1] = (n+1)*sk;
        SFToV[sf++][2] = FACE_OPENBS;
    }
    /* right - outflow */
    for(n=0;n<My;n++){
        SFToV[sf  ][0] = Mx+n*sk;
        SFToV[sf  ][1] = Mx+(n+1)*sk;
        SFToV[sf++][2] = FACE_OPENBS;
    }
#if DEBUG
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    FILE *fp = create_log(__FUNCTION__, nprocs, procid);
    print_int_matrix2file(fp,"SFToV", SFToV, Nsurf, 3);
#endif

    dg_grid_add_BS2d(grid, Nsurf, SFToV);
#if DEBUG
    print_int_matrix2file(fp, "EToV", grid->EToV, grid->K, dg_cell_Nfaces(grid->cell));
    print_int_matrix2file(fp, "EToBS", grid->EToBS, grid->K, dg_cell_Nfaces(grid->cell));
    fclose(fp);
#endif
    matrix_int_free(SFToV);
    return;
}

static dg_phys* rotation_phys_init(dg_grid *grid){
    dg_region *region = dg_region_create(grid);
    dg_mesh *mesh = dg_mesh_create(region);
    dg_edge *edge = dg_edge_create(mesh);
    dg_phys *phys = dg_phys_create(3, edge);
    const int K = dg_grid_K(dg_phys_grid(phys));
    const int Np = dg_cell_Np( dg_phys_cell(phys));
    const double sigma = 125*1e3/(33*33);
    const double w = 5*M_PI/6;
    const double xc = 0.0;
    const double yc = 0.6;

    double **x = dg_region_x(dg_phys_region(phys));
    double **y = dg_region_y(dg_phys_region(phys));

    dg_real *f_Q = dg_phys_f_Q(phys);

    int k,n,sk = 0;
    for(k=0;k<K;++k){
        for(n=0;n<Np;++n){
            const double xt = x[k][n];
            const double yt = y[k][n];
            double t = -sigma * (( xt - xc )*( xt - xc ) + ( yt - yc )*( yt - yc ));
            f_Q[sk++] = (dg_real) exp(t); // c field
            f_Q[sk++] = (dg_real)(-w * yt); // flow rate at x-coordinate
            f_Q[sk++] = (dg_real)( w * xt); // flow rate at y-coordinate
        }
    }
    return phys;
}

static dg_phys* advection_phys_init(dg_grid *grid){
    dg_region *region = dg_region_create(grid);
    dg_mesh *mesh = dg_mesh_create(region);
    dg_edge *edge = dg_edge_create(mesh);
    dg_phys *phys = dg_phys_create(3, edge);

    const int K = dg_grid_K(dg_phys_grid(phys));
    const int Np = dg_cell_Np(dg_phys_cell(phys));
    const int Nfield = dg_phys_Nfield(phys);
    const double xc = -0.5;
    const double yc = -0.5;
    //const double sigma = 125*1e3/(33*33);

    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    double **x = dg_region_x(dg_phys_region(phys));
    double **y = dg_region_y(dg_phys_region(phys));
    extern Conv_Solver solver;
    // read input file
    arg_section **sec_p = conv_read_inputfile(solver.filename);
    /// 4. phys info
    arg_section *sec = sec_p[4];
    double ut,vt,miu;
    sscanf(sec->arg_vec_p[0], "%lf\n", &ut);
    sscanf(sec->arg_vec_p[1], "%lf\n", &vt);
    sscanf(sec->arg_vec_p[2], "%lf\n", &miu);
    if(!procid){
        printf(HEADLINE " u: %lf\n", ut);
        printf(HEADLINE " v: %lf\n", vt);
        printf(HEADLINE " miu: %lf\n", miu);
        if(miu < EPS){
            fprintf(stderr, "The viscosity value %lf should be positive.\n", miu);
            exit(-1);
        }
    }
    if(miu > 1e-6) { solver.vis_flag = 1; }

    conv_arg_section_free(sec_p);
    dg_real *f_Q = dg_phys_f_Q(phys);

    int k,n;
    for(k=0;k<K;++k){
        for(n=0;n<Np;++n){
            int sk = (k*Np+n)*Nfield;
            const double xt = x[k][n];
            const double yt = y[k][n];
            double t = -( ( xt - xc )*( xt - xc ) + ( yt - yc )*( yt - yc ) )/miu;
            f_Q[sk+0] = (dg_real) exp(t); // c field
            f_Q[sk+1] = (dg_real) ut; // flow rate at x-coordinate
            f_Q[sk+2] = (dg_real) vt; // flow rate at y-coordinate
        }
    }

    dg_phys_LDG *ldg = phys->ldg;
    for(k=0;k<K;k++){
        for(n=0;n<Np;n++){
            int sk = (k*Np+n)*Nfield;
            dg_phys_ldg_sqrt_miux(ldg)[sk] = miu;
            dg_phys_ldg_sqrt_miuy(ldg)[sk] = miu;
        }
    }

#if DEBUG
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    FILE *fp = create_log(__FUNCTION__, nprocs, procid);
    print_double_vector2file(fp, "miux", dg_phys_ldg_sqrt_miux(ldg), K*Np*Nfield);
    print_double_vector2file(fp, "miuy", dg_phys_ldg_sqrt_miuy(ldg), K*Np*Nfield);
    fclose(fp);
#endif
    return phys;
}

static dg_phys* pure_advection_phys_init(dg_grid *grid){

    dg_region *region = dg_region_create(grid);
    dg_mesh *mesh = dg_mesh_create(region);
    dg_edge *edge = dg_edge_create(mesh);
    dg_phys *phys = dg_phys_create(3, edge);

    const int K = dg_grid_K(dg_phys_grid(phys));
    const int Np = dg_cell_Np(dg_phys_cell(phys));
    const double xc = -0.5;
    const double yc = -0.5;
    const double sigma = 125*1e3/(33*33);
    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    double **x = dg_region_x(dg_phys_region(phys));
    double **y = dg_region_y(dg_phys_region(phys));
    extern Conv_Solver solver;
    // read input file
    arg_section **sec_p = conv_read_inputfile(solver.filename);
    /// 4. phys info
    arg_section *sec = sec_p[4];
    double ut,vt;
    sscanf(sec->arg_vec_p[0], "%lf\n", &ut);
    sscanf(sec->arg_vec_p[1], "%lf\n", &vt);
    if(!procid){
        printf(HEADLINE " u: %lf\n", ut);
        printf(HEADLINE " v: %lf\n", vt);
    }
    solver.vis_flag = 0;
    conv_arg_section_free(sec_p);
    dg_real *f_Q = dg_phys_f_Q(phys);

    int k,n,sk = 0;
    for(k=0;k<K;++k){
        for(n=0;n<Np;++n){
            const double xt = x[k][n];
            const double yt = y[k][n];
            double t = -sigma * (( xt - xc )*( xt - xc ) + ( yt - yc )*( yt - yc ));
            f_Q[sk++] = (dg_real) exp(t); // c field
            f_Q[sk++] = (dg_real)ut; // flow rate at x-coordinate
            f_Q[sk++] = (dg_real)vt; // flow rate at y-coordinate
        }
    }
    return phys;
}

static void set_time_info(double dt){
    extern Conv_Solver solver;
    arg_section **sec_p = conv_read_inputfile(solver.filename);

    double cfl,dt_user,ftime,out_dt;
    arg_section *sec = sec_p[3];
    sscanf(sec->arg_vec_p[0], "%lf\n", &(cfl));
    sscanf(sec->arg_vec_p[1], "%lf\n", &(dt_user));
    sscanf(sec->arg_vec_p[2], "%lf\n", &(ftime));

    /// 5. output time interval
    sec = sec_p[5];
    sscanf(sec->arg_vec_p[1], "%lf\n", &(out_dt));
    conv_arg_section_free(sec_p);

    //if(dt_user < dt) {dt = dt_user;}
    dt = dt_user;
    double gdt;
    MPI_Allreduce(&dt, &gdt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt = cfl*gdt;

    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    if(!procid) printf(HEADLINE " cfl: %f\n", cfl);
    if(!procid) printf(HEADLINE " dt: %f\n", dt);
    if(!procid) printf(HEADLINE " final time: %f\n", ftime);
    if(!procid) printf(HEADLINE " output dt: %f\n", out_dt);

    solver.finaltime = ftime;
    solver.dt = dt;
    solver.outDt = out_dt;
    return;
}

static void advdiff_time_init(dg_phys *phys){
    double dt = 1e6; // initial time step
    const int N = dg_cell_N(dg_phys_cell(phys));
    const int Np = dg_cell_Np(dg_phys_cell(phys));
    const int K = dg_grid_K(dg_phys_grid(phys));
    double *len = dg_phys_region(phys)->len;

    dg_phys_LDG *ldg = phys->ldg;
    const double miu = dg_phys_ldg_sqrt_miux(ldg)[0];

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
            dt = min(dt, r*r/sqrt(miu));
        }
    }

    set_time_info(dt);
    return;
}

static void advection_time_init(dg_phys *phys){
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
    set_time_info(dt);
    return;
}

