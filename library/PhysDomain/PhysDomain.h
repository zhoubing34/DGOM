//
// Created by li12242 on 16/6/20.
//

#ifndef DGOM_PHYSDOMAIN_H
#define DGOM_PHYSDOMAIN_H

#include "MultiRegions/MultiRegions.h"

typedef struct PhysDomain2d{
    /** number of variable fields */
    int Nfields;
    /** mesh */
    MultiReg2d* mesh;
    /** surf infomation */
    float *surfinfo;
    /** variable array */
    float *f_Q, *f_rhsQ, *f_resQ;

    /** send/recv data */
    float *f_inQ, *f_outQ;
}PhysDomain2d;

#endif //DGOM_PHYSDOMAIN_H
