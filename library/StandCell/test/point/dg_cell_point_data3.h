//
// Created by li12242 on 17/3/15.
//

#ifndef DGOM_DG_CELL_LINE_DATA3_H
#define DGOM_DG_CELL_LINE_DATA3_H

#define NP 1
#define NFP 0
#define NV 1

double point_r[NP] = {0};

double point_V[NP][NP] = {{1}};

double point_M[NP][NP] = {{1}};

double point_Dr[NP][NP] = {{ 0 }};

double point_LIFT[NP][NFP];

double point_VX[NV] = {0};

#endif //DGOM_DG_CELL_LINE_DATA3_H
