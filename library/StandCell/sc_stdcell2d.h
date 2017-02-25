//
// Created by li12242 on 17/2/25.
//

#ifndef DGOM_SC_STDCELL2D_H
#define DGOM_SC_STDCELL2D_H

#include "sc_stdcell.h"

/* get the gradient matrix Dr and Ds of Lagrange basis at (r,s) at order N */
void sc_deriMatrix2d(stdCell *cell, void (*derorthfunc)(stdCell *, int ind, double *dr, double *ds));
/* calculate the Gauss quadrature weights for faces (ws) and volume (wv) integral */
void sc_GaussQuadrature2d(stdCell *cell);
/* create mass matrix of edges */
void sc_surfMassMatrix2d(stdCell *cell, double **Mes);

#endif //DGOM_SC_STDCELL2D_H
