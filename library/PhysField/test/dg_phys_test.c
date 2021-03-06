//
// Created by li12242 on 12/11/16.
//
#include "dg_phys_test.h"
#include "dg_phys_strong_vol_opt_test.h"
#include "dg_phys_strong_surf_opt_test.h"
#include "dg_phys_obc_test.h"
#include "dg_phys_limiter_test.h"
#include "dg_phys_strong_LDG_opt_test.h"

int Mx = 4;
int My = 4;
int N = 3;
int Nfield = 2;
double xmin = -1, xmax = 1;
double ymin = -1, ymax = 1;

/** generation of physfield */
dg_phys *uniform_tri_physfield();
dg_phys *uniform_quad_physfield();
dg_phys *rectangle_tri_physfield();

#define Nphys 3
typedef dg_phys* (*phys_func)();
static const phys_func phys_creator[Nphys] = {
        uniform_tri_physfield,
        uniform_quad_physfield,
        rectangle_tri_physfield,
};

/**
 * @brief create uniform triangle mesh of physField
 * @return physField
 */
dg_phys *uniform_tri_physfield(){
    int type = 1;
    dg_cell *tri = dg_cell_creat(N, TRIANGLE);
    dg_area *area = dg_area_create_uniform(tri, Mx, My, xmin, xmax, ymin, ymax, type);
    dg_phys *phys = dg_phys_create(Nfield, area);
    return phys;
}
/**
 * @brief create uniform quadrilateral mesh of physField
 * @return physField
 */
dg_phys *uniform_quad_physfield(){
    dg_cell *quad = dg_cell_creat(N, QUADRIL);
    int type = 1;
    dg_area *area = dg_area_create_uniform(quad, Mx, My, xmin, xmax, ymin, ymax, type);
    dg_phys *phys = dg_phys_create(Nfield, area);
    return phys;
}

dg_phys *rectangle_tri_physfield(){
    dg_cell *tri = dg_cell_creat(N, TRIANGLE);
    char casename[] = "SWE2d/Rectangle/Rectangle";
    dg_area *area = dg_area_create_from_file(tri, casename);
    dg_phys *phys = dg_phys_create(2, area);
    char obcfile[] = "SWE2d/Rectangle/Rectangle.obc4.nc";
    phys->attach_obc_ncfile(phys, obcfile);
    return phys;
}

/**
 * @brief deallocate the memory of physField
 * @param phys
 */
void phys_free(dg_phys *phys){
    dg_cell_free(dg_phys_cell(phys));
    dg_grid_free(dg_phys_grid(phys));
    dg_region_free(dg_phys_region(phys));
    dg_mesh_free(dg_phys_mesh(phys));
    dg_edge_free(dg_phys_edge(phys));
    dg_phys_free(phys);
}

/** test functions */
#define Ntest 5
typedef int (*test_func)(dg_phys *, int verbose);
static const test_func phys_test_func[Ntest] = {
        dg_phys_strong_vol_opt2d_test,
        dg_phys_strong_surf_opt2d_test,
        dg_phys_obc_test,
        dg_phys_limiter_test,
        dg_phys_strong_LDG_opt2d_test,
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
    test_command(argc, argv, &ishelp, &isverbose);

    // print help information
    if(ishelp){
        char helpstr[] = HEADEND "DGOM\n" HEADLINE "Unit tests for PhysField library\n"
                HEADLINE "Example usages:\n"
                HEADLINE "   mpirun -n 2 -host localhost ./physfield_test\n"
                HEADLINE "\n"
                HEADLINE "Optional features:\n"
                HEADLINE "   -help     print help information\n"
                HEADEND  "   -verbose  print variables to log files\n\n";
        printf("%s",helpstr);
        return 0;
    }

    MPI_Init(&argc, &argv);
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if(!procid){
        printf(HEADSTART "Running %d test for %d dg_phys from %s.\n",
               Ntest, Nphys, __FUNCTION__);
    }

    dg_phys *phys;
    int n,m,tmp;

    for(n=0;n<Nphys;n++){
        /* create dg_phys */
        phys = phys_creator[n]();
        if(!procid){ printf(HEADLINE "%d test for dg_phys[%d]\n", Ntest, n); }

        int Nfail=0;
        for(m=0;m<Ntest;m++){
            tmp = phys_test_func[m](phys, isverbose);
            if(tmp){ Nfail++; }
        }
        phys_free(phys);
        if(Nfail){
            if(!procid){ printf(HEADEND "%d test failed for dg_phys[%d]\n\n", Nfail, n); }
        } else{
            if(!procid){ printf(HEADEND "%d test passed for dg_phys[%d]\n\n", Ntest, n);}
        }
    }


    MPI_Finalize();

    return 0;
}