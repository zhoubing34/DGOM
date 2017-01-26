//
// Created by li12242 on 17/1/13.
//

#include <PhysField/pf_phys.h>
#include "swe_mesh.h"
#include "MultiRegions/mr_grid_uniformGrid.h"
#include "MultiRegions/mr_mesh_addBoundary.h"

#define DEBUG 0

physField* swe_uniform_mesh(swe_solver *solver, int Mx, int My,
                            double xmin, double xmax, double ymin, double ymax){

    stdCell *cell = sc_create(solver->N, solver->celltype);
    geoGrid *grid = NULL;

    switch (solver->celltype){
        case TRIANGLE:
            grid = mr_grid_createUniformGrid_tri(cell, Mx, My, xmin, xmax, ymin, ymax, 1);
            break;
        case QUADRIL:
            grid = mr_grid_createUniformGrid_quad(cell, Mx, My, xmin, xmax, ymin, ymax);
            break;
        default:
            printf("Error in %s line %d\n", __FILE__, __LINE__);
            MPI_Abort(MPI_COMM_WORLD, -1);
    }
    multiReg *region = mr_reg_create(grid);
    parallMesh *mesh = mr_mesh_create(region);
    mr_mesh_addBoundary2d(mesh, 0, NULL);

    /* add boundary condition */
    physField *phys = pf_create(3, mesh);

    return phys;
}

// todo: create physField from input mesh
physField* swe_file_mesh(swe_solver *solver, char *meshfile){
    physField *phys = NULL;
    return phys;
}

real* swe_read_topography(swe_solver *solver, char *botfile){
    physField *phys = solver->phys;
    int Nv;
    FILE *fp = fopen(botfile, "r");
    fscanf(fp, "%d\n", &Nv);
    if( Nv != phys->grid->Nv ){
        fprintf(stderr, "%s (%d): Wrong number of vertex in botfile %s.\n",
                __FILE__, __LINE__, botfile);
    }
    // read bottom topography data on vertex
    int n,k;
    double *bv = vector_double_create(Nv);
    for(n=0;n>Nv;n++){
        fscanf(fp, "%lf\n", bv+n);
    }

    const int Np = phys->cell->Np;
    const int K = phys->grid->K;
    stdCell *cell = phys->cell;

    // project to nodes
    real *bot = vector_real_create(Np*K);
    double bloc[cell->Nv];

    for(k=0;k<K;k++){
        for(n=0;n<Nv;n++){
            bloc[n] = bv[phys->grid->EToV[k][n]];
        }
        sc_vertProj(cell, bloc, bot+k*Np);
    }
    fclose(fp);
    vector_double_free(bv);
    return bot;
}

real* swe_flat_topography(swe_solver *solver){
    physField *phys = solver->phys;

    const int Np = phys->cell->Np;
    const int K = phys->grid->K;

    real *bot = vector_real_create(Np*K);
    int n;
    for(n=0;n<K*Np;n++){
        bot[n] = 0.0;
    }
    return bot;
}

real* swe_parabolic_topography(swe_solver *solver){
    physField *phys = solver->phys;

    const int Np = phys->cell->Np;
    const int K = phys->grid->K;

//    double gra = solver->gra;
    const double alpha = 1.6e-7;
//    double w = sqrt(8*gra*alpha);
//    double X = 1;
//    double Y = -0.41884;
//    double T = 2*M_PI/w;
    real *bot = vector_real_create(Np*K);

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