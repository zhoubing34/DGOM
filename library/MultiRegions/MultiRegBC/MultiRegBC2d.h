//
// Created by li12242 on 16/12/14.
//

#ifndef DGOM_MULTIREGBC2D_H
#define DGOM_MULTIREGBC2D_H

#include "MultiRegions/MultiRegions.h"

typedef enum {
    INNERLOC, // 0-inner face
    INNERBS,  // 1-inner boundary surface
    OPENBS,   // 2-open boundary
    SLIPWALL, // 3-slip wall
    NSLIPWALL // 4-non-slip wall
} BCType;

typedef struct {
    int **EToBS; // element to boundary surface type
    int Nobs;    // # of open boundary surface
    real *vert_ext; // sparse vector for external data on vertex

    /* send & recv data */
    real *f_ext; // array of external data for open boundary
    real *f_in;  //
    real *f_out; //

}MultiRegBC2d;

#endif //DGOM_MULTIREGBC2D_H
