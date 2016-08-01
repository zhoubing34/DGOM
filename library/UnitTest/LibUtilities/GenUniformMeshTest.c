#include "LibUtiltiesTest.h"

int main(int argc, char **argv){
    double xmin = -1, xmax = 1;
    double ymin = -1, ymax = 1;

    int **EToV;
    double *VX, *VY;
    int Nv, Ne;

    /* tri */
    GenUniformTriMesh(1, 1, xmin, xmax, ymin, ymax, 0, &Ne, &Nv, &EToV, &VX, &VY);
    PrintIntMatrix("EToV for triangle", EToV, Ne, 3);
    PrintVector("VX for triangle", VX, Nv);
    PrintVector("VY for triangle", VY, Nv);

    DestroyIntMatrix(EToV);
    DestroyVector(VX);
    DestroyVector(VY);

    /* quad */
    GenUniformQuadMesh(2, 2, xmin, xmax, ymin, ymax, &Ne, &Nv, &EToV, &VX, &VY);
    PrintIntMatrix("EToV for quad", EToV, Ne, 4);
    PrintVector("VX for quad", VX, Nv);
    PrintVector("VY for quad", VY, Nv);

    DestroyIntMatrix(EToV);
    DestroyVector(VX);
    DestroyVector(VY);

    return 0;
}

