#ifndef DGOM_SETTESTMULTIREGIONS_H
#define DGOM_SETTESTMULTIREGIONS_H

#include "MultiRegions/MultiRegions.h"

MultiReg2d* SetTriMultiRegions();
MultiReg2d* SetQuadMultiRegions();
MultiReg2d* SetTriParallelMultiRegions();
MultiReg2d* SetQuadParallelMultiRegions();

#endif //DGOM_SETTESTMULTIREGIONS_H