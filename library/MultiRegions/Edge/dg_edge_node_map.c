//
// Created by li12242 on 17/3/10.
//

#include "dg_edge_node_map.h"

void dg_edge_node_map3d(dg_edge *edge){
    return;
}

void dg_edge_node_map2d(dg_edge *edge){

    dg_cell *cell = edge->cell;
    const int Nedge = dg_edge_Nedge(edge);
    const int Np = dg_cell_Np(cell);
    const int TotalParNode = dg_mesh_Nparn(edge->mesh);
    int **Fmask = cell->Fmask;
    double **x = edge->region->x;
    double **y = edge->region->y;
    int *parnode = edge->mesh->parnode;
    // count total nodes on edges
    int f,totalNode = 0;
    for(f=0;f<Nedge;f++){
        int f1 = edge->varfM[f];
        totalNode += dg_cell_Nfp(cell, f1);
    }
    // map node
    int *varpM = vector_int_create(totalNode);
    int *varpP = vector_int_create(totalNode);

    int n1,n2,sk=0;
    for(f=0;f<Nedge;f++){
        int k1 = edge->varkM[f];
        int k2 = edge->varkP[f];
        int f1 = edge->varfM[f];
        int f2 = edge->varfP[f];
        int ftype = edge->ftype[f];

        const int Nfp = dg_cell_Nfp(cell, f1);
        for(n1=0;n1<Nfp;n1++){
            const int idM = k1*Np + Fmask[f1][n1];
            varpM[sk] = idM;
            double xM = x[0][idM];
            double yM = y[0][idM];

            if(ftype != INNERBS){
                for(n2=0;n2<Nfp;n2++){
                    const int idP = k2*Np + Fmask[f2][n2];
                    double xP = x[0][idP];
                    double yP = y[0][idP];
                    double d12 = (xM-xP)*(xM-xP) + (yM-yP)*(yM-yP);
                    if(d12 < EPS){ varpP[sk] = idP; break; }
                }
            }else{ // parallel face
                for(n2=0;n2<TotalParNode;n2++){
                    if( varpM[sk] == parnode[n2] ){ varpP[sk] = n2; break; }
                }
            }
            sk++;
        }
    }
    // assignment
    edge->Nnode = totalNode;
    edge->varpM = varpM;
    edge->varpP = varpP;

    return;
}