//
// Created by li12242 on 17/1/13.
//

#include <PhysField/pf_phys.h>
#include "swe_mesh.h"
#include "MultiRegions/Grid/dg_grid_reader.h"
#include "MultiRegions/Mesh/mr_mesh_bc.h"

#define DEBUG 0
#if DEBUG
#include "Utility/UTest.h"
#endif

physField* swe_uniform_mesh(swe_solver *solver, int Mx, int My,
                            double xmin, double xmax, double ymin, double ymax){

    dg_cell *cell = dg_cell_creat(solver->N, solver->celltype);
    dg_grid *grid = NULL;

    switch (solver->celltype){
        case TRIANGLE:
            grid = dg_grid_create_uniform_tri(cell, Mx, My, xmin, xmax, ymin, ymax, 1);
            break;
        case QUADRIL:
            grid = dg_grid_create_uniform_quad(cell, Mx, My, xmin, xmax, ymin, ymax);
            break;
        default:
            printf("Error in %s line %d\n", __FILE__, __LINE__);
            MPI_Abort(MPI_COMM_WORLD, -1);
    }
    dg_region *region = mr_reg_create(grid);
    dg_mesh *mesh = mr_mesh_create(region);
    mr_mesh_add_bc2d(mesh, 0, NULL);

    /* add boundary condition */
    physField *phys = pf_create(3, mesh);

    return phys;
}

static void swe_boundary_condition(swe_solver *solver, dg_mesh *mesh){
    mr_mesh_read_bcfile2d(mesh, solver->casename);
}

/**
 * @brief create physical field from input mesh file
 * */
physField* swe_file_mesh(swe_solver *solver, char *meshfile){

    dg_cell *cell = dg_cell_creat(solver->N, solver->celltype);
    dg_grid *grid = dg_grid_read_file2d(cell, meshfile);
    dg_region *region = mr_reg_create(grid);
    dg_mesh *mesh = mr_mesh_create(region);
    swe_boundary_condition(solver, mesh);
    physField *phys = pf_create(3, mesh);
#if DEBUG
    char casename[20] = "swe_triGrid_test";
    FILE *fp = CreateLog(casename, grid->procid, grid->nprocs);
    PrintIntMatrix2File(fp, "EToV", grid->EToV, grid->K, grid->cell->Nv);
    PrintVector2File(fp, "vx", grid->vx, grid->Nv);
    PrintVector2File(fp, "vy", grid->vy, grid->Nv);
    PrintMatrix2File(fp, "J", region->J, grid->K, grid->cell->Np);
    fclose(fp);
#endif
    return phys;
}

dg_real* swe_read_topography(swe_solver *solver, char *botfile){
    physField *phys = solver->phys;
    int Nvert;
    FILE *fp;
    if( (fp = fopen(botfile, "r")) == NULL ){
        fprintf(stderr, "swe_read_topography (%s): %d\n"
                        "Unable to open topography file %s.\n",
                __FILE__,__LINE__,botfile);
    }
    fscanf(fp, "%d\n", &Nvert);
    if( Nvert != phys->grid->Nv ){
        fprintf(stderr, "%s (%d): Wrong number of vertex in botfile %s.\n",
                __FILE__, __LINE__, botfile);
    }
    // read bottom topography data on vertex
    int n,k,tmp;
    double *bv = vector_double_create(Nvert);
    for(n=0;n<Nvert;n++){
        fscanf(fp, "%d", &tmp);
        fscanf(fp, "%lf", bv+n);
    }
    const int Np = phys->cell->Np;
    const int K = phys->grid->K;
    dg_cell *cell = phys->cell;

    // project to nodes
    dg_real *bot = vector_real_create(Np*K);
    const int Nv = cell->Nv;
    double bloc[Nv];

    for(k=0;k<K;k++){
        for(n=0;n<Nv;n++){
            bloc[n] = bv[phys->grid->EToV[k][n]];
        }
        dg_cell_proj_vert2node(cell, bloc, bot + k * Np);
    }
    fclose(fp);
    vector_double_free(bv);
    return bot;
}

dg_real* swe_flat_topography(swe_solver *solver){
    physField *phys = solver->phys;

    const int Np = phys->cell->Np;
    const int K = phys->grid->K;

    dg_real *bot = vector_real_create(Np*K);
    int n;
    for(n=0;n<K*Np;n++){
        bot[n] = 0.0;
    }
    return bot;
}

dg_real* swe_parabolic_topography(swe_solver *solver){
    physField *phys = solver->phys;

    const int Np = phys->cell->Np;
    const int K = phys->grid->K;

//    double gra = solver->gra;
    const double alpha = 1.6e-7;
//    double w = sqrt(8*gra*alpha);
//    double X = 1;
//    double Y = -0.41884;
//    double T = 2*M_PI/w;
    dg_real *bot = vector_real_create(Np*K);

    int n,k,sk=0;
    for(k=0;k<K;k++){
        for(n=0;n<Np;n++){
            double x = phys->region->x[k][n];
            double y = phys->region->y[k][n];
            bot[sk++] = alpha*(x*x+y*y);
        }
    }
    return bot;
}