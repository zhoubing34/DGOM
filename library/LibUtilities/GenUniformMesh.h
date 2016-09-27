#ifndef UNIFORMMESH_H
#define UNIFORMMESH_H

#include "LibUtilities.h"
#include "LibUtilities/UnstructMesh.h"

UnstructMesh* GenUniformTriMesh(
        int Mx, int My,
        double xmin, double xmax,
        double ymin, double ymax, int type);
UnstructMesh* GenUniformQuadMesh(
        int Mx, int My,
        double xmin, double xmax,
        double ymin, double ymax);
UnstructMesh* GenParallelUniformTriMesh(
        int Mx, int My,
        double xmin, double xmax,
        double ymin, double ymax, int type,
        int procid, int nprocs);
UnstructMesh* GenParallelUniformQuadMesh(
        int Mx, int My,
        double xmin, double xmax,
        double ymin, double ymax,
        int procid, int nprocs);

#endif // UNIFORMMESH_H
