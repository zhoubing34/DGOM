//
// Created by li12242 on 16/12/14.
//

#ifndef DGOM_MULTIREGBC2D_H
#define DGOM_MULTIREGBC2D_H

#include "MultiRegions/MultiRegions.h"

typedef enum {
    INNERLOC, // 0-inner face
    INNERBS,  // 1-inner boundary surface
    SLIPWALL, // 2-slip wall
    NSLIPWALL, // 3-non-slip wall
    OPENBS   // 4-open boundary
} BCType;


typedef struct{
    int Nv;
    int *BVToV;
} vertlist;

typedef struct {
    MultiReg2d *mesh;
    /* adjacent cell id of each cell */
    int **EToE;
    /* adjacent face type of each cell */
    int **EToBS;
    /* array of boundary type indicators */
    int *bcTypeList;

    /** volume id of -ve trace of face node */
    int *vmapM;
    /** volume id of +ve trace of face node */
    int *vmapP;

    /** total number of nodes to send recv */
    int parNodeTotalOut;
    /** index list of nodes to send out */
    int *nodeIndexOut;

    /* total number of element's value to send and recv */
    int parCellTotalOut;
    /* index list of elements to send out */
    int *cellIndexOut;

    int Nobc;    // number of open boundary surface
    vertlist **oblist; // vector of open boundary pointer

}MultiRegBC2d;

MultiRegBC2d* MultiRegBC2d_create(MultiReg2d *mesh, int Nsurf, int **SFToV);
void MultiRegBC2d_free(MultiRegBC2d* bc2d);

#endif //DGOM_MULTIREGBC2D_H
