//
// Created by li12242 on 16/12/27.
//

#include "conv_driver2d.h"
#include "conv_mesh.h"
#include "MultiRegions/mr_grid_create.h"
#include "MultiRegions/mr_mesh_bc.h"

/* set open boundary condition */
static void conv_setOBC(int *indicator, int **SFToV);

/* set parallel mesh object */
parallMesh* conv_mesh(stdCell *shape){

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
            printf("Error in %s line %d\n", __FILE__, __LINE__);
            MPI_Abort(MPI_COMM_WORLD, -1);
    }

    multiReg *region = mr_reg_create(grid);
    parallMesh *mesh = mr_mesh_create(region);

    /* add boundary condition */
    int indicator[4] = {5,5,5,5};
    int Nsurf = 2*(Ne+Ne);
    int **SFToV = matrix_int_create(Nsurf, 3);
    conv_setOBC(indicator, SFToV);

    mr_mesh_addBoundary2d(mesh, Nsurf, SFToV);

    matrix_int_free(SFToV);

    return mesh;
}


static void conv_setOBC(int *indicator, int **SFToV){

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