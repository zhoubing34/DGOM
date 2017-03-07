#ifndef MR_MESH_PRARLLELPAIRS
#define MR_MESH_PRARLLELPAIRS

/* ParallelPairs.c */
void mr_mesh_parallelPairs(void *objs, int Nmyobjs, int sizeobj,
                           int  (*numget)(const void *),
                           void (*numset)(const void *, int),
                           int  (*procget)(const void *),
                           void (*marry)(const void *, const void *),
                           int (*compare_objs)(const void *, const void *));

#endif