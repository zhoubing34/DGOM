#ifndef DGOM_PHYSDOMAIN_H
#define DGOM_PHYSDOMAIN_H

#include "MultiRegions/mr_mesh.h"

typedef struct{
    int dim; ///< dimensions
    int Nfield; ///< number of variable fields

    parallMesh *mesh; ///< parallel mesh object
    multiReg *region; // multi-region object
    geoGrid *grid; ///< geometry grid
    stdCell *cell; ///< standard element

    int Nsurfinfo; ///< number of elements in surfinfo
    real *surfinfo; ///< surface infomation

    int Nvgeo; ///< number of elements in vgeo
    real *vgeo; ///< volume geometry information

    /* parallel buffer and map */
    int parallNodeNum; ///< number of node variables to send/recv
    int *nodeIndexOut; ///< map from send array `f_outQ` to node index */
    real *f_inQ, *f_outQ; ///< send/recv buffers for nodal information

    int parallCellNum; ///< number of cell variables to send/recv
    int *cellIndexOut; ///< map from send buffer `c_outQ` to cell index
    int *cellIndexIn; ///< map from recv buffer `c_inQ` to cell index
    real *c_inQ, *c_outQ; ///< send/recv buffers for elemental information

    /* information */
    real *c_Q; ///< elemental information
    real *f_ext; ///< external data

    real *f_Q, *f_rhsQ, *f_resQ; ///< nodal information
} physField;

physField* phys_create(int Nfields, parallMesh *mesh);
void PhysDomain2d_free(physField *phys);

#endif //DGOM_PHYSDOMAIN_H