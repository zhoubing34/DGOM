//
// Created by li12242 on 17/1/29.
//

#ifndef DGOM_MAT_ARITH_H
#define DGOM_MAT_ARITH_H

#ifdef MKL_LIB
#include "mkl_lapacke.h"
#else
#include "f2c.h"
#include "blaswrap.h"
#include "clapack.h"
#endif

#endif //DGOM_MAT_ARITH_H
