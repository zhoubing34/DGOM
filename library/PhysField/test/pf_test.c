//
// Created by li12242 on 12/11/16.
//


#include "MultiRegions/mr_grid_create.h"
#include "MultiRegions/mr_mesh_bc.h"
#include "pf_test.h"
#include "pf_nodeFetch_test.h"
#include "pf_strong_volume_flux2d_test.h"
#include "pf_strong_surface_flux2d_test.h"
#include "pf_strong_viscosity_LDG_flux2d_test.h"
#include "pf_cellMean_test.h"
#include "pf_cellFetch_test.h"
#include "pf_limiter_test.h"

int Mx = 4;
int My = 2;
int N = 3;
int Nfield = 2;
double xmin = -1, xmax = 1;
double ymin = -1, ymax = 1;

/** generation of physfield */
physField *uniform_tri_nobc_physfield();
physField *uniform_quad_nobc_physfield();

#define Nphys 2
typedef physField* (*phys_func)();
static const phys_func phys_creator[Nphys] = {
        uniform_tri_nobc_physfield,
        uniform_quad_nobc_physfield
};

/**
 * @brief create uniform triangle mesh of physField
 * @return physField
 */
physField *uniform_tri_nobc_physfield(){
    int type = 1;
    stdCell *tri = sc_create(N, TRIANGLE);
    geoGrid *tri_grid = mr_grid_create_uniform_tri(tri, Mx, My, xmin, xmax, ymin, ymax, type);
    multiReg *tri_region = mr_reg_create(tri_grid);
    parallMesh *tri_mesh = mr_mesh_create(tri_region);
    mr_mesh_addBoundary2d(tri_mesh, 0, NULL);
    physField *tri_phys = pf_create(Nfield, tri_mesh);
    return tri_phys;
}
/**
 * @brief create uniform quadrilateral mesh of physField
 * @return physField
 */
physField *uniform_quad_nobc_physfield(){
    stdCell *quad = sc_create(N, QUADRIL);
    geoGrid *quad_grid = mr_grid_create_uniform_quad(quad, Mx, My, xmin, xmax, ymin, ymax);
    multiReg *quad_region = mr_reg_create(quad_grid);
    parallMesh *quad_mesh = mr_mesh_create(quad_region);
    mr_mesh_addBoundary2d(quad_mesh, 0, NULL);
    physField *quad_phys = pf_create(Nfield, quad_mesh);
    return quad_phys;
}

/**
 * @brief deallocate the memory of physField
 * @param phys
 */
void phys_nobc_free(physField *phys){
    /* free memory */
    sc_free(phys->cell);
    mr_grid_free(phys->grid);
    mr_reg_free(phys->region);
    mr_mesh_free(phys->mesh);
    pf_free(phys);
}

/** test functions */
#define Ntest 7
typedef int (*test_func)(physField *, int verbose);
static const test_func phys_test_func[Ntest] = {
        phys_cellMean_test,
        phys_nodeFetch_test,
        phys_cellFetch_test,
        phys_strong_volume_flux2d_test,
        phys_strong_surface_flux2d_test,
        phys_strong_viscosity_LDG_flux2d_test,
        phys_limiter_test
};

/**
 * @brief test for physField model
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv){
    // get settings
    int ishelp, isverbose;
    UTest_Command(argc, argv, &ishelp, &isverbose);

    // print help information
    if(ishelp){
        printf(HEADEND "DGOM\n" HEADLINE "Unit tests for PhysField library\n"
                       HEADLINE "Example usages:\n"
                       HEADLINE "   mpirun -n 2 -host localhost ./physfield_test\n"
                       HEADLINE "\n"
                       HEADLINE "Optional features:\n"
                       HEADLINE "   -help     print help information\n"
                       HEADEND  "   -verbose  print variables to log files\n\n");
        return 0;
    }

    MPI_Init(&argc, &argv);
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if(!procid){
        printf(HEADSTART "Running %d test for %d physField from %s.\n",
               Ntest, Nphys, __FILE__);
    }

    physField *phys;
    int n,m,tmp;

    for(n=0;n<Nphys;n++){
        /* create physField */
        phys = phys_creator[n]();
        if(!procid){ printf(HEADLINE "%d test for physField[%d]\n", Ntest, n); }

        int failNum=0;
        for(m=0;m<Ntest;m++){
            tmp = phys_test_func[m](phys, isverbose);
            if(tmp){ failNum++; }
        }
        phys_nobc_free(phys);
        if(failNum){
            if(!procid){ printf(HEADEND "%d test failed for physField[%d]\n\n", failNum, n); }
        } else{
            if(!procid){ printf(HEADEND "%d test passed for physField[%d]\n\n", Ntest, n);}
        }
    }


    MPI_Finalize();

    return 0;
}