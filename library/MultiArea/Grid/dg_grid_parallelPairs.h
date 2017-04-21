//
// Created by li12242 on 17/3/8.
//

#ifndef DGOM_DG_GRID_PARALLELPAIRS_H
#define DGOM_DG_GRID_PARALLELPAIRS_H

/* ParallelPairs.c */
void dg_grid_parallelPairs(void *objs, int Nmyobjs, int sizeobj,
                           int  (*numget)(const void *),
                           void (*numset)(const void *, int),
                           int  (*procget)(const void *),
                           void (*marry)(const void *, const void *),
                           int (*compare_objs)(const void *, const void *));


#endif //DGOM_DG_GRID_PARALLELPAIRS_H
