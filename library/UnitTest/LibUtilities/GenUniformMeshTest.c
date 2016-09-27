#include "LibUtiltiesTest.h"

int main(int argc, char **argv){
    double xmin = -1, xmax = 1;
    double ymin = -1, ymax = 1;

    UnstructMesh *grid;

    /* tri */
    grid = UniformTriMesh_create(1, 1, xmin, xmax, ymin, ymax, 0);
    PrintIntMatrix("EToV for triangle", grid->EToV, grid->ne, 3);
    PrintVector("VX for triangle", grid->vx, grid->nv);
    PrintVector("VY for triangle", grid->vy, grid->nv);

    UnstructMesh_free(grid);

    /* quad */
    grid = UniformQuadMesh_create(2, 2, xmin, xmax, ymin, ymax);
    PrintIntMatrix("EToV for quad", grid->EToV, grid->ne, 4);
    PrintVector("VX for quad", grid->vx, grid->nv);
    PrintVector("VY for quad", grid->vy, grid->nv);

    UnstructMesh_free(grid);

    return 0;
}

