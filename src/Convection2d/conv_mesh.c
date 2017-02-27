//
// Created by li12242 on 16/12/27.
//

#include "conv_driver2d.h"
#include "conv_mesh.h"
#include "MultiRegions/mr_grid_create.h"
#include "MultiRegions/mr_mesh_bc.h"

/* set open boundary condition */
static void conv_uniform_bc(int *indicator, int **SFToV);
static parallMesh* conv_uniform_mesh(stdCell *shape);
static parallMesh* conv_userset_mesh(stdCell *shape);

/* set parallel mesh object */
parallMesh* conv_mesh(stdCell *shape){

    parallMesh *mesh = NULL;

    extern conv_solver2d solver;
    switch (solver.caseid){
        case conv_advection_diffusion:
            mesh = conv_uniform_mesh(shape); break;
        case conv_rotational_convection:
            mesh = conv_uniform_mesh(shape); break;
        case conv_userset:
            mesh = conv_userset_mesh(shape); break;
        default:
            fprintf(stderr, "%s: %d\nUnknown test case %d\n",
                    __FUNCTION__, __LINE__, solver.caseid);
            exit(-1);
    }
    return mesh;
}

static parallMesh* conv_userset_mesh(stdCell *shape){
    extern conv_solver2d solver;
    geoGrid *grid = mr_grid_read_file2d(shape, solver.casename);
    multiReg *region = mr_reg_create(grid);
    parallMesh *mesh = mr_mesh_create(region);
    mr_mesh_read_bcfile2d(mesh, solver.casename);
    return mesh;
}

static parallMesh* conv_uniform_mesh(stdCell *shape){
    extern conv_solver2d solver;
    const int Ne = solver.Ne;
    geoGrid *grid = NULL;
    switch (solver.celltype){
        case TRIANGLE:
            grid = mr_grid_create_uniform_tri(shape, Ne, Ne, -1, 1, -1, 1, 1);
            break;
        case QUADRIL:
            grid = mr_grid_create_uniform_quad(shape, Ne, Ne, -1, 1, -1, 1);
            break;
        default:
            fprintf(stderr, "%s: %d\nUnknown cell type %d\n",
                    __FUNCTION__, __LINE__, solver.celltype);
            MPI_Abort(MPI_COMM_WORLD, -1);
    }
    multiReg *region = mr_reg_create(grid);
    parallMesh *mesh = mr_mesh_create(region);

    // boundary condition
    /* add boundary condition */
    int indicator[4] = {5,5,5,5};
    int Nsurf = 2*(Ne+Ne);
    int **SFToV = matrix_int_create(Nsurf, 3);
    conv_uniform_bc(indicator, SFToV);
    mr_mesh_add_bc2d(mesh, Nsurf, SFToV);
    matrix_int_free(SFToV);
    return mesh;
}

static void conv_uniform_bc(int *indicator, int **SFToV){

    extern conv_solver2d solver;
    const int Mx=solver.Ne, My=solver.Ne;

    int Nfaces[4] = {Mx, Mx, My, My};
    int start_face_Id[4] = {0, Mx, Mx+Mx, Mx*2+My};
    int start_vert_Id[4] = {0, (Mx+1)*My, Mx, 0};
    int vert_stride[4] = {1, 1, Mx+1, Mx+1};
    // loop over 4 boundaries
    int b,f;
    for(b=0;b<4;b++){
        int ind = indicator[b];
        int startid = start_vert_Id[b];
        int vstride = vert_stride[b];
        int startf = start_face_Id[b];
        for(f=0;f<Nfaces[b];f++){
            SFToV[f+startf][0] = startid; startid+=vstride;
            SFToV[f+startf][1] = startid;
            SFToV[f+startf][2] = ind;
        }
    }
    return;
}