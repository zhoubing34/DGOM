#include "UnstructMesh.h"

/**
 * @brief
 * Deallocate unstructed grid structure.
 */
void DestroyUnstructMesh(UnstructMesh *grid) {
    free(grid->EToV[0]);
    free(grid->EToV);
    switch (grid->dim) {
        case 1:
            DestroyVector(grid->vx);
            break;
        case 2:
            DestroyVector(grid->vx);
            DestroyVector(grid->vy);
            break;
        case 3:
            DestroyVector(grid->vx);
            DestroyVector(grid->vy);
            DestroyVector(grid->vz);
            break;
        default:
            printf("Wrong dimension %d of unstructed grid %s.", grid->dim, grid->name);
            exit(-1);
    }
    return;
}

/**
 * @brief
 * Create and allocate an unstructed grid of single element type.
 */
UnstructMesh* CreateUnstructMesh(int Dim, int Ne, int Nv, EleType eletype){

    UnstructMesh *grid = (UnstructMesh*) calloc(1, sizeof(UnstructMesh));

    grid->dim = Dim;
    grid->ne  = Ne;
    grid->nv  = Nv;

    int perNv;
    // the vertex num in each element
    switch (eletype) {
        case TRIANGLE:
            perNv=3; break;
        case QUADRIL:
            perNv=4; break;
        case TETRA:
            perNv=4; break;
        case TRIPRISM:
            perNv=6; break;
        case HEXA:
            perNv=8; break;
        default:
            printf("Unknown element type!");
            exit(-1);
    }
    // allocate
    grid->EToV=BuildIntMatrix(Ne, perNv);
    switch (Dim){ // allocate vertex
        case 1:
            grid->vx  =BuildVector(Nv);
            break;
        case 2:
            grid->vx  =BuildVector(Nv);
            grid->vy  =BuildVector(Nv);
            break;
        case 3:
            grid->vx  =BuildVector(Nv);
            grid->vy  =BuildVector(Nv);
            grid->vz  =BuildVector(Nv);
            break;
        default:
            printf("Wrong dimension %d of unstructed grid.", Dim);
            exit(-1);
    }
    return grid;
}