//
// Created by li12242 on 12/19/16.
//

#include "multiregion_test.h"
#include "dg_reg_test.h"
#include "dg_mesh_test.h"
#include "dg_grid_test.h"
#include "dg_edge_test.h"

dg_grid *uniform_tri_grid();
dg_grid *uniform_quad_grid();
dg_grid *user_set_tri_grid();

/** global variable */
int Mx = 4;
int My = 4;
int N  = 2;
double xmin = -1, xmax = 1;
double ymin = -1, ymax = 1;

/** mesh creating function */
#define Nmesh 1
typedef dg_grid* (* grid_create_func)();
grid_create_func grid_func[Nmesh] = {uniform_quad_grid};

/**
 * @brief create uniform triangle mesh
 * @return
 */
dg_grid *uniform_tri_grid(){
    dg_cell *shape = dg_cell_creat(N, TRIANGLE);
    dg_grid *grid = dg_grid_uniform_tri(shape, Mx, My, xmin, xmax, ymin, ymax, 1);
    return grid;
}
/**
 * @brief create uniform quadrilateral mesh
 * @return
 */
dg_grid *uniform_quad_grid(){
    dg_cell *shape = dg_cell_creat(N, QUADRIL);
    dg_grid *grid = dg_grid_uniform_quad(shape, Mx, My, xmin, xmax, ymin, ymax);
    return grid;
}
/**
 * @brief create triangle mesh from files
 * @return
 */
dg_grid *user_set_tri_grid(){
    dg_cell *shape = dg_cell_creat(N, TRIANGLE);
    char casename[] = "example/SWE2d/Rectangle/Rectangle";
    dg_grid *grid = dg_grid_read_file2d(shape, casename);
    dg_grid_read_BSfile2d(grid, casename);
    return grid;
}

/**
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv){

    // get settings
    int ishelp, isverbose;
    test_command(argc, argv, &ishelp, &isverbose);

    MPI_Init(&argc, &argv);
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // print help information
    if(ishelp && !procid){
        printf(HEADEND "DGOM\n" HEADLINE "Unit tests for MultiRegions library\n"
                       HEADLINE "Usages:\n"
                       HEADLINE "   mpirun -n 2 -host localhost ./MulitRegions_Test\n"
                       HEADLINE "\n"
                       HEADLINE "Optional features:\n"
                       HEADLINE "   -help     print help information\n"
                       HEADEND  "   -verbose  print variables to log files\n");
        return 0;
    }

    // local variable
    int err, n;
    // test
    dg_mesh *mesh;
    dg_region *region;
    dg_grid *grid;
    dg_edge *edge;
    for(n=0;n<Nmesh;n++){
        if(!procid){ printf(HEADSTART "Running test for mesh[%d]\n", n); }

        int Nfail = 0;
        //mesh = mesh_func[n]();
        grid = grid_func[n]();

        err = dg_grid_Ktol_test(grid, isverbose); if(err){ Nfail+=1; }
        err = dg_grid_vertex_test(grid, isverbose); if(err){ Nfail+=1; }
        err = dg_grid_EToV_test(grid, isverbose); if(err){ Nfail+=1; }
        err = dg_grid_connect_test(grid, isverbose); if(err){ Nfail+=1; }
        err = dg_grid_EToBS_test(grid, isverbose); if(err){ Nfail+=1; }

        region = dg_region_create(grid);
        err = dg_region_node_test(region, isverbose); if(err){ Nfail+=1; }
        err = dg_region_volume_factor_test(region, isverbose); if(err){ Nfail+=1; }
        err = dg_region_face_factor_test(region, isverbose); if(err){ Nfail+=1; }
        err = dg_region_scale_test(region, isverbose); if(err){ Nfail+=1; }
        err = dg_region_face_integral_test(region, isverbose); if(err){ Nfail+=1; }

        mesh = dg_mesh_create(region);
        err = dg_mesh_parallel_test(mesh, isverbose); if(err){ Nfail+=1; }
        err = dg_mesh_cell_fetch_buffer_test(mesh, isverbose); if(err){ Nfail+=1; }
        err = dg_mesh_node_fetch_buffer_test(mesh, isverbose); if(err){ Nfail+=1; }

        edge = dg_edge_create(mesh);
        err = dg_edge_facemap_test(edge, isverbose); if(err){ Nfail+=1; }
        err = dg_edge_nodemap_test(edge, isverbose); if(err){ Nfail+=1; }
        err = dg_edge_surfinfo_test(edge, isverbose); if(err){ Nfail+=1; }

        dg_grid_free(grid);
        dg_region_free(region);
        dg_mesh_free(mesh);
        dg_edge_free(edge);

        if(Nfail) {
            if(!procid) printf(HEADEND "%d test faild for mesh[%d]\n", Nfail, n);
        } else {
            if(!procid) printf(HEADEND "all test passed for mesh[%d]\n", n);
        }
    }

    MPI_Finalize();

    return 0;
}

