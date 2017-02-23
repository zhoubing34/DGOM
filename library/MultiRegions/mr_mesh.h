//
// Created by li12242 on 12/19/16.
//

#ifndef DGOM_MR_MESH_H
#define DGOM_MR_MESH_H

#include "mr_reg.h"

typedef enum {
    INNERLOC=0, // 0-inner face, adjacent face in the same process
    INNERBS=1,  // 1-inner boundary surface, adjacent face in other process
    SLIPWALL=2, // 2-slip wall
    NSLIPWALL=3, // 3-non-slip wall
    OPENBS   // 4-open boundary
} mr_surfType;

typedef struct {

    int dim; ///< dimension
    int procid; ///< process id
    int nprocs; ///< number of process

    multiReg *region; ///< multi-region object
    geoGrid *grid; ///< geometry grid object
    stdCell *cell; ///< standard element object

    /* adjacent cell id of each cell */
    int **EToE; ///< adjacent element id in each element
    int **EToP; ///< adjacent element process
    int **EToF; ///< adjacent face index (0,1,2...)

    int *Npar; ///< number of faces adjacent to each process

    int *vmapM; ///< node id of -ve trace of face node
    int *vmapP; ///< node id of +ve trace of face node

    int parallNodeNum; ///< total number of nodes to send recv
    int *nodeIndexOut; ///< index list of nodes to send out

    int parallCellNum; ///< total number of elements to send and recv
    int *cellIndexOut; ///< index list of elements in send buffer
    int *cellIndexIn; ///< index list of elements in recv buffer
    int *faceIndexIn; ///< index list of face in recv buffer

    /* boudary condition */
    int **EToBS; ///< adjacent face type of each cell
    int Nbc; ///< number of boundary surface (except INNERLOC and INNERBS)
    int *bcind; ///< array of boundary type indicators
    int Nobc;    ///< number of open boundary surface
    int *obcind; ///< array of open boundary type indicators
    char **obcfilename; ///< filename of open boundary files

} parallMesh;

parallMesh* mr_mesh_create(multiReg *region);
void mr_mesh_free(parallMesh *mesh);


#endif //DGOM_MR_MESH_H
