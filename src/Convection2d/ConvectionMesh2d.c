#include "ConvectionDriver2d.h"
#include "LibUtilities/GenUniformMesh.h"

/* private function */
MultiReg2d* ReadTriMesh(StdRegions2d *shape, int Ne);
MultiReg2d* ReadQuadMesh(StdRegions2d *shape, int Ne);


MultiReg2d* ReadMesh(StdRegions2d *shape, int Ne){
    MultiReg2d *mesh;

    if (shape->Nfaces == 3){
        mesh = ReadTriMesh(shape, Ne);
    }else if(shape->Nfaces == 4){
        mesh = ReadQuadMesh(shape, Ne);
    }else{
        fprintf(stderr, "Wrong number of element shapes, Nfaces=%d \n", shape->Nfaces);
        exit(-1);
    }
    return mesh;
}

/**
 * @brief
 * Generate local uniform triangle mesh for each processor
 *
 * @details
 * The whole domain is [-1, 1]x[-1, 1], with each edge partition to *Ne* 
 * elements. The vertex list is arranged from upper left to lower right
 * directions. The whole mesh is generate with newEToV for element to 
 * vertex list. Each process reads part of mesh into mesh->EToV and its 
 * corresponding vertex coordinate mesh->GX and mesh->GY. The uniform 
 * triangle element is separated from the square, and the separation is 
 * controlled by UP_RIGHT macro.
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 * @return
 * return values：
 * name     | type      | description of value
 * -------- |---------- |----------------------
 * mesh     | Mesh*     | mesh object
 *
 */
MultiReg2d* ReadTriMesh(StdRegions2d *shape, int Ne){

    /* decide on parition */
    int procid, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* generate uniform grid */
    int type=0; // 0 for '\' and 1 for '/'
    UnstructMesh *triGrid = ParallelUniformTriMesh_create(
            Ne, Ne, -1, 1, -1, 1, type, procid, nprocs);

    /* MultiRegions object and allocation */
    MultiReg2d *region = MultiReg2d_create(shape, triGrid);

    /* deallocate grid */
    UnstructMesh_free(triGrid);

    return region;
}

/**
 * @brief
 * Generate uniform quadrilateral mesh
 *
 * @details
 * The whole domain is [-1, 1]x[-1, 1], with each edge partition to *Ne* 
 * elements. The vertex list is arranged from upper left to lower right
 * directions. The whole mesh is generate with newEToV for element to 
 * vertex list. Each processor reads part of mesh into mesh->EToV and 
 * its corresponding vertex coordinate mesh->GX and mesh->GY.
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 * @return
 * return values：
 * name     | type      | description of value
 * -------- |---------- |----------------------
 * mesh     | Mesh*     | mesh object
 *
 */
MultiReg2d* ReadQuadMesh(StdRegions2d *shape, int Ne){

    /* decide on parition */
    int procid, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* generate uniform grid */
    UnstructMesh *quadGrid = ParallelUniformQuadMesh_create(
            Ne, Ne, -1, 1, -1, 1, procid, nprocs);

    /* MultiRegions object and allocation */
    MultiReg2d *region = MultiReg2d_create(shape, quadGrid);

    /* deallocate grid */
    UnstructMesh_free(quadGrid);

    return region;
}

