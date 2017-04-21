//
// Created by li12242 on 16/12/27.
//
#include "conv2d.h"

#define DEBUG 0

static dg_area* user_grid_init();
static dg_phys* user_phys_init(dg_area *area);
static void conv_time_init(dg_phys *phys);
static void set_const_vis_value(dg_phys *phys, double vis);

void conv_init(){
    dg_area *area = user_grid_init();
    dg_phys *phys = user_phys_init(area);
    conv_time_init(phys);
    return;
}

static dg_area* user_grid_init(){
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
    dg_area *area = dg_area_create_from_file(cell, casename);

    return area;
}

static dg_phys* user_phys_init(dg_area *area){
    extern Conv_Solver solver;
    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    dg_phys *phys = dg_phys_create(3, area);
    /// 2. obc file
    extern Conv_Solver solver;
    arg_section **sec_p = conv_read_inputfile(solver.filename);
    char parameter_str[MAX_NAME_LENGTH];
    arg_section *sec = sec_p[2];
    strcpy(parameter_str, sec->arg_vec_p[0]);
    if(!procid) printf(HEADLINE " open boundary file: %s\n", parameter_str);
    if(strlen(parameter_str)) { phys->attach_obc_ncfile(phys, parameter_str); }

    /// 4. initial condition
    sec = sec_p[4];
    strcpy(parameter_str, sec->arg_vec_p[0]);
    if(!procid) printf(HEADLINE " initial condition file: %s\n", parameter_str);
    phys->initialize_from_file(phys, parameter_str);

    strcpy(parameter_str, sec->arg_vec_p[1]);
    solver.vis_flag = 0;
    if(strlen(parameter_str)) {
        if (!procid) printf(HEADLINE " viscosity value: %s\n", parameter_str);
        double vis;
        sscanf(sec->arg_vec_p[1], "%lf\n", &(vis));
        set_const_vis_value(phys, vis);
        solver.vis_flag = 1;
    }

    conv_arg_section_free(sec_p);
    /* assignment */
    solver.phys = phys;

    return phys;
}

static void set_const_vis_value(dg_phys *phys, double vis){
    const int K = dg_grid_K(dg_phys_grid(phys));
    const int Nfield = dg_phys_Nfield(phys);
    const int Np = dg_cell_Np(dg_phys_cell(phys));

    dg_phys_LDG *ldg = phys->ldg;
    const double sqrt_miu = sqrt(vis);
    int k,n,fld;
    for(k=0;k<K;k++){
        for(n=0;n<Np;n++){
            for(fld=0;fld<Nfield;fld++){
                int sk = (k*Np+n)*Nfield+fld;
                dg_phys_ldg_sqrt_miux(ldg)[sk] = sqrt_miu;
                dg_phys_ldg_sqrt_miuy(ldg)[sk] = sqrt_miu;
            }
        }
    }
    return;
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

    double cfl,dt_user,ftime,out_dt;
    arg_section *sec = sec_p[3];
    sscanf(sec->arg_vec_p[0], "%lf\n", &(cfl));
    sscanf(sec->arg_vec_p[1], "%lf\n", &(dt_user));
    sscanf(sec->arg_vec_p[2], "%lf\n", &(ftime));

    /// 5. output time interval
    sec = sec_p[5];
    sscanf(sec->arg_vec_p[0], "%lf\n", &(out_dt));

    conv_arg_section_free(sec_p);

    if(dt_user < dt) {dt = dt_user;}
    double gdt;
    MPI_Allreduce(&dt, &gdt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt = cfl*gdt;

    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    if(!procid) printf(HEADLINE " dt: %f\n", dt);
    if(!procid) printf(HEADLINE " final time: %f\n", ftime);
    if(!procid) printf(HEADLINE " output dt: %f\n", out_dt);
    // assignment
    extern Conv_Solver solver;
    solver.finaltime = ftime;
    solver.dt = dt;
    solver.outDt = out_dt;
    return;
}
