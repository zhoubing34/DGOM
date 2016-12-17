#ifndef UNIFORMMESH_H
#define UNIFORMMESH_H

#include "MultiRegions/UnstructMesh.h"

UnstructMesh* UniformTriMesh_create(
        int Mx, int My,
        double xmin, double xmax,
        double ymin, double ymax, int type);
UnstructMesh* UniformQuadMesh_create(
        int Mx, int My,
        double xmin, double xmax,
        double ymin, double ymax);
UnstructMesh* ParallelUniformTriMesh_create(
        int Mx, int My,
        double xmin, double xmax,
        double ymin, double ymax, int type,
        int procid, int nprocs);
UnstructMesh* ParallelUniformQuadMesh_create(
        int Mx, int My,
        double xmin, double xmax,
        double ymin, double ymax,
        int procid, int nprocs);

#endif // UNIFORMMESH_H
