//
// Created by li12242 on 12/19/16.
//

#include <MultiRegions/mr_mesh.h>
#include "mr_test.h"
#include "mr_reg_test.h"
#include "mr_mesh_test.h"
#include "mr_grid_test.h"

parallMesh *uniform_tri_mesh();
parallMesh *uniform_quad_mesh();
parallMesh *user_set_tri_mesh();

/** global variable */
int Mx = 4;
int My = 2;
int N  = 2;
double xmin = -1, xmax = 1;
double ymin = -1, ymax = 1;

/** mesh creating function */
#define Nmesh 3
typedef parallMesh* (* mesh_create_func)();
mesh_create_func mesh_func[Nmesh] = {uniform_tri_mesh, uniform_quad_mesh, user_set_tri_mesh};

/**
 * @brief create uniform triangle mesh
 * @return
 */
parallMesh *uniform_tri_mesh(){
    dg_cell *shape = dg_cell_creat(N, TRIANGLE);
    dg_grid *grid = mr_grid_create_uniform_tri(shape, Mx, My, xmin, xmax, ymin, ymax, 1);
    multiReg *region = mr_reg_create(grid);
    parallMesh *mesh = mr_mesh_create(region);
    mr_mesh_add_bc2d(mesh, 0, NULL);
    return mesh;
}
/**
 * @brief create uniform quadrilateral mesh
 * @return
 */
parallMesh *uniform_quad_mesh(){
    dg_cell *shape = dg_cell_creat(N, QUADRIL);
    dg_grid *grid = mr_grid_create_uniform_quad(shape, Mx, My, xmin, xmax, ymin, ymax);
    multiReg *region = mr_reg_create(grid);
    parallMesh *mesh = mr_mesh_create(region);
    mr_mesh_add_bc2d(mesh, 0, NULL);
    return mesh;
}
/**
 * @brief create triangle mesh from files
 * @return
 */
parallMesh *user_set_tri_mesh(){
    dg_cell *shape = dg_cell_creat(N, TRIANGLE);
    char casename[] = "example/Rectangle/tri/Rectangle";
    dg_grid *grid = mr_grid_read_file2d(shape, casename);
    multiReg *region = mr_reg_create(grid);
    parallMesh *mesh = mr_mesh_create(region);
    mr_mesh_read_bcfile2d(mesh, casename);
    return mesh;
}

void mesh_free(parallMesh *mesh){
    dg_cell_free(mesh->cell);
    mr_grid_free(mesh->grid);
    mr_mesh_del_bc2d(mesh);
    mr_mesh_free(mesh);
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
    parallMesh *mesh;
    multiReg *region;
    dg_grid *grid;
    for(n=0;n<Nmesh;n++){
        if(!procid){ printf(HEADSTART "Running test for mesh[%d]\n", n); }

        int Nfail = 0;
        mesh = mesh_func[n]();
        region = mesh->region;
        grid = mesh->grid;

        err = mr_grid_Ktol_test(grid, isverbose); if(err){ Nfail+=1; }
        err = mr_grid_vertex_test(grid, isverbose); if(err){ Nfail+=1; }
        err = mr_grid_EToV_test(grid, isverbose); if(err){ Nfail+=1; }

        err = mr_reg_node_test(region, isverbose); if(err){ Nfail+=1; }
        err = mr_reg_volume_factor_test(region, isverbose); if(err){ Nfail+=1; }
        err = mr_reg_face_factor_test(region, isverbose); if(err){ Nfail+=1; }
        err = mr_reg_scale_test(region, isverbose); if(err){ Nfail+=1; }

        err = mr_mesh_connet_test(mesh, isverbose); if(err){ Nfail+=1; }
        err = mr_mesh_parallel_test(mesh, isverbose); if(err){ Nfail+=1; }
        err = mr_mesh_bc_test(mesh, isverbose); if(err){ Nfail+=1; }

        mesh_free(mesh);

        if(Nfail) {
            if(!procid) printf(HEADEND "%d test faild for mesh[%d]\n", Nfail, n);
        } else {
            if(!procid) printf(HEADEND "all test passed for mesh[%d]\n", n);
        }
    }

    MPI_Finalize();

    return 0;
}

