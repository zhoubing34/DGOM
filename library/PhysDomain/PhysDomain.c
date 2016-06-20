/**
 * @file
 * Set physical variables
 *
 * @brief
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "PhysDomain.h"

PhysDomain2d* GetPhysDomain2d(MultiReg2d *mesh, int Nfields){
    PhysDomain2d *phys = (PhysDomain2d *) calloc(1, sizeof(PhysDomain2d));

    /* number of variables */
    phys->Nfields = Nfields;

    /* mesh */
    phys->mesh = mesh;
    /* data array */
    phys->f_Q    = (float *) calloc(mesh->K*mesh->StdElement->Np*Nfields, sizeof(float));
    phys->f_rhsQ = (float *) calloc(mesh->K*mesh->StdElement->Np*Nfields, sizeof(float));
    phys->f_resQ = (float *) calloc(mesh->K*mesh->StdElement->Np*Nfields, sizeof(float));

    /* MPI send/recv buffer */
    phys->f_inQ  = (float *) calloc(mesh->parNtotalout*Nfields, sizeof(float));
    phys->f_outQ = (float *) calloc(mesh->parNtotalout*Nfields, sizeof(float));

    phys->surfinfo = (float *) calloc(mesh->K*)

}
