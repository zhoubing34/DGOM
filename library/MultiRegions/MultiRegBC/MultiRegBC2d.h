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
} OBC2d;

typedef struct {
    int **EToBS; // element to boundary surface type
    int Nobc;    // number of open boundary surface
    int *bcTypeList;
    real *vert_ext; // sparse vector for external data on vertex

    OBC2d **obc2d; // vector of open boundary pointer

//    /* send & recv data */
//    real *f_ext; // array of external data for open boundary
//    real *f_in;  //
//    real *f_out; //

}MultiRegBC2d;

#endif //DGOM_MULTIREGBC2D_H
