#include "UnstructMesh.h"

/**
 * @brief
 * Deallocate unstructured grid structure.
 */
void UnstructMesh_free(UnstructMesh *grid) {
    free(grid->EToV[0]);
    free(grid->EToV);
    switch (grid->dim) {
        case 1:
            Vector_free(grid->vx);
            break;
        case 2:
            Vector_free(grid->vx);
            Vector_free(grid->vy);
            break;
        case 3:
            Vector_free(grid->vx);
            Vector_free(grid->vy);
            Vector_free(grid->vz);
            break;
        default:
            printf("Wrong dimension %d of unstructed grid %s.", grid->dim, grid->name);
            exit(-1);
    }
    return;
}

/**
 * @brief
 * Create and allocate an unstructured grid of single element type.
 */
UnstructMesh* UnstructMesh_create(int Dim, int Ne, int Nv, ElementType eletype){

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
    grid->EToV= IntMatrix_create(Ne, perNv);
    switch (Dim){ // allocate vertex
        case 1:
            grid->vx  = Vector_create(Nv);
            break;
        case 2:
            grid->vx  = Vector_create(Nv);
            grid->vy  = Vector_create(Nv);
            break;
        case 3:
            grid->vx  = Vector_create(Nv);
            grid->vy  = Vector_create(Nv);
            grid->vz  = Vector_create(Nv);
            break;
        default:
            printf("Wrong dimension %d of unstructed grid.", Dim);
            exit(-1);
    }
    return grid;
}