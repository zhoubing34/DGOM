#ifndef DGOM_PHYSDOMAIN_H
#define DGOM_PHYSDOMAIN_H

#include "MultiRegions/mr_mesh.h"


typedef struct {
    real *px_Q; ///< dfdx partial derivative for x direction
    real *py_Q; ///< dfdy partial derivative for y direction
    real *px_inQ, *px_outQ; ///< send and recv buffers for p_Q
    real *py_inQ, *py_outQ; ///< send and recv buffers for q_Q
    real *vis_Q; ///< viscosity on each node
} pf_LDG_solver;


typedef struct{
    int dim; ///< dimensions
    int Nfield; ///< number of variable fields

    parallMesh *mesh; ///< parallel mesh object
    multiReg *region; ///< multi-region object
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
    //int *cellIndexIn; ///< map from recv buffer `c_inQ` to cell index
    real *c_inQ, *c_outQ; ///< send/recv buffers for elemental information

    /* information */
    real *c_Q; ///< elemental information
    real *f_ext; ///< external data

    real *f_Q; ///< nodal information
    real *f_rhsQ; ///< RHS data
    real *f_resQ; ///< residual data

    pf_LDG_solver *viscosity; ///< LDG solver
} physField;

physField* pf_create(int Nfields, parallMesh *mesh);
void pf_free(physField *phys);

#endif //DGOM_PHYSDOMAIN_H
