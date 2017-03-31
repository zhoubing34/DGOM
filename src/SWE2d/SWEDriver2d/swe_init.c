#include "swe2d.h"
#include "../SWELib/swe_lib.h"

#define DEBUG 0

static dg_grid *swe_user_init_grid();
static dg_phys* user_phys_init(dg_grid *grid);
static void swe_read_manning_frict_file(dg_phys *phys, char *filename);

typedef struct SWE_Init_Creator{
    dg_grid *(*grid_init)();
    dg_phys *(*phys_init)(dg_grid *grid);
}SWE_Init_Creator;

const static SWE_Init_Creator user_init = {swe_user_init_grid, user_phys_init};

void swe_init(){
    const SWE_Init_Creator *creator = &user_init;

    extern SWE_Solver solver;
    dg_grid *grid = creator->grid_init();
    dg_phys *phys = creator->phys_init(grid);

    solver.phys = phys;
    return;
}

static dg_grid *swe_user_init_grid(){
    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    arg_section **sec_p = swe_read_section();
    /// 0. case info
    char casename[MAX_NAME_LENGTH];
    arg_section *sec = sec_p[0];
    strcpy(casename, sec->arg_vec_p[0]);

    if(!procid){ printf(HEAD_LINE " casename: %s\n", casename); }
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
                printf(HEAD_LINE " cell type: triangle\n"); break;
            case QUADRIL:
                printf(HEAD_LINE " cell type: quadrilateral\n"); break;
            default:
                fprintf(stderr, "%s (%d)\nUnknown cell type %d.\n",
                        __FUNCTION__, __LINE__, cell_type);
                MPI_Abort(MPI_COMM_WORLD, -1);
        }
        printf(HEAD_LINE " polynomial degree: %d\n", N);
    }
    swe_free_section(sec_p);
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
    dg_phys *phys = dg_phys_create(4, edge); // physical field for (h, hu, hv, b)
    /// 2. obc file
    arg_section **sec_p = swe_read_section();
    char obcfile[MAX_NAME_LENGTH];
    arg_section *sec = sec_p[2];
    strcpy(obcfile, sec->arg_vec_p[0]);
    if(strlen(obcfile)) { phys->obc_add(phys, obcfile); }
    if(!procid) printf(HEAD_LINE " open boundary file: %s\n", obcfile);

    /// 3. time info
    extern SWE_Solver solver;
    sec = sec_p[3];
    sscanf(sec->arg_vec_p[0], "%lf\n", &(solver.cfl));
    sscanf(sec->arg_vec_p[1], "%lf\n", &(solver.dt));
    sscanf(sec->arg_vec_p[2], "%lf\n", &(solver.ftime));

    if(!procid) printf(HEAD_LINE " cfl: %lf\n", solver.cfl);
    if(!procid) printf(HEAD_LINE " dt: %lf\n", solver.dt);
    if(!procid) printf(HEAD_LINE " final time: %lf\n", solver.ftime);


    /// 4. initial condition
    char filename[MAX_NAME_LENGTH];
    sec = sec_p[4];
    sscanf(sec->arg_vec_p[0], "%lf\n", &(solver.gra)); /* 0. read gra */
    sscanf(sec->arg_vec_p[1], "%lf\n", &(solver.hcrit)); /* 1. read hcrit */
    strcpy(filename, sec->arg_vec_p[2]); /* 2. read Manning friction coefficient */
    swe_read_manning_frict_file(phys, filename);
    strcpy(filename, sec->arg_vec_p[3]); /* 3. read initial condition file */
    phys->init_file(phys, filename);

    if(!procid) printf(HEAD_LINE " gra: %lf\n", solver.gra);
    if(!procid) printf(HEAD_LINE " hcrit: %lf\n", solver.hcrit);
    if(!procid) printf(HEAD_LINE " Manning friction file: %s\n", filename);
    if(!procid) printf(HEAD_LINE " initial condition file: %s\n", filename);

    /// 5. output info
    sec = sec_p[5];
    strcpy(solver.outfilename, sec->arg_vec_p[0]); /* 0. output file name */
    sscanf(sec->arg_vec_p[1], "%lf\n", &(solver.outDt));

    if(!procid) printf(HEAD_LINE " output file name: %s\n", solver.outfilename);
    if(!procid) printf(HEAD_LINE " output time interval: %lf\n", solver.outDt);

    swe_free_section(sec_p);
    return phys;
}

static void swe_read_manning_frict_file(dg_phys *phys, char *filename){
    FILE *fp;
    dg_fopen(fp, filename, "Unable to open manning coefficient file");

    dg_cell *cell = dg_phys_cell(phys);
    const int K = dg_grid_K(dg_phys_grid(phys));
    const int Np = dg_cell_Np(cell);
    const int Nvert = dg_grid_Nv(dg_phys_grid(phys)); ///< number of vertex in grid;
    const int Nv = dg_cell_Nv(cell); ///< number of vertex in cell;
    int k,n,nvert;
    fscanf(fp, "%d\n", &nvert); // read vertex number
    if(Nvert != nvert ){
        fprintf(stderr, "%s (%d): The vertex node number in %s is incorrect\n",
                __FUNCTION__, __LINE__, filename);
        exit(-1);
    }
    dg_real *vert_M = vector_real_create(Nvert);
    for(n=0;n<Nvert;n++){
        int tmp; // node index
        fscanf(fp, "%d %lf", &tmp, vert_M+n);
    }

    /* read element file */
    int **EToV = dg_grid_EToV(dg_phys_grid(phys));
    extern SWE_Solver solver;
    dg_real *Mann = vector_real_create(Np*K);

    dg_real node_Mann[Np], vert_Mann[Nv];
    for(k=0;k<K;k++) {
        // vertex initial values
        for (n=0;n<Nv;n++) { vert_Mann[n] = vert_M[EToV[k][n]]; }
        cell->proj_vert2node(cell, 1, vert_Mann, node_Mann);
        // assignment
        for(n=0;n<Np;n++){ Mann[k*Np+n] = node_Mann[n]; }
    }

    solver.m = Mann;
    vector_real_free(vert_M);
    fclose(fp);
    return;
}

//static void swe_dambreakdry_init(SWE_Solver *solver){
//
//    physField *phys = solver->phys;
//    const int K = phys->grid->K;
//    const int Np = phys->cell->Np;
//
//    const double x0 = 500; // dam position
//    register int k,n,sk=0;
//
//    dg_real *f_Q = phys->f_Q;
//    for(k=0;k<K;k++){
//        double *xt = phys->region->x[k];
//        double area = phys->region->size[k];
//        double xc = mr_reg_integral(phys->region, k, xt)/area;
//
//        for(n=0;n<Np;n++){
//            // h field
//            if(xc<x0){
//                f_Q[sk++] = 10.0;
//            }else{
//                f_Q[sk++] = 0.0;
//            }
//            f_Q[sk++] = 0.0; // q_x field
//            f_Q[sk++] = 0.0; // q_y field
//        }
//    }
//    return;
//}
//
//static void swe_dambreakwet_init(SWE_Solver *solver){
//
//    physField *phys = solver->phys;
//    const int K = phys->grid->K;
//    const int Np = phys->cell->Np;
//
//    const double x0 = 500; // dam position
//    register int k,n,sk=0;
//
//    dg_real *f_Q = phys->f_Q;
//    for(k=0;k<K;k++){
//        double *xt = phys->region->x[k];
//        double area = phys->region->size[k];
//        double xc = mr_reg_integral(phys->region, k, xt)/area;
//
//        for(n=0;n<Np;n++){
//            // h field
//            if(xc<x0){
//                f_Q[sk++] = 10.0;
//            }else{
//                f_Q[sk++] = 2.0;
//            }
//            f_Q[sk++] = 0.0; // q_x field
//            f_Q[sk++] = 0.0; // q_y field
//        }
//    }
//    return;
//}
//
//static void swe_parabolicbowl_init(SWE_Solver *solver){
//
//    physField *phys = solver->phys;
//    const int K = phys->grid->K;
//    const int Np = phys->cell->Np;
//
//    const double alpha = 1.6e-7;
//    const double X = 1;
//    const double Y = -0.41884;
//
//    dg_real *f_Q = phys->f_Q;
//    register int k,n,sk=0;
//    for(k=0;k<K;k++){
//        for(n=0;n<Np;n++){
//            const double xt = phys->region->x[k][n];
//            const double yt = phys->region->y[k][n];
//            const double r2 = (xt*xt + yt*yt);
//            const double r2ext = (X+Y)/(alpha*(X*X - Y*Y));
//            if(r2<r2ext){
//                f_Q[sk++] = 1.0/(X+Y) + alpha*(Y*Y - X*X)*r2/(X+Y)/(X+Y);
//                f_Q[sk++] = 0.0;
//                f_Q[sk++] = 0.0;
//            }else{
//                f_Q[sk++] = 0.0;
//                f_Q[sk++] = 0.0;
//                f_Q[sk++] = 0.0;
//            }
//        }
//    }
//    return;
//}