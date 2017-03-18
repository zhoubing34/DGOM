//
// Created by li12242 on 17/3/10.
//

#include "dg_edge_node_map.h"

void dg_edge_node_map3d(dg_edge *edge){
    return;
}

void dg_edge_node_map(dg_edge *edge){

    dg_cell *cell = edge->cell;
    dg_region *region = edge->region;
    const int Nedge = dg_edge_Nedge(edge);
    const int Np = dg_cell_Np(cell);
    const int TotalParNode = dg_mesh_NfetchNode(edge->mesh);
    const int Nfaces = dg_cell_Nfaces(cell);
    int **Fmask = dg_cell_Fmask(cell);
    double **x = dg_region_x(region);
    double **y = dg_region_y(region);
    double **z = dg_region_z(region);

    int *parnode = edge->mesh->NBFToN;
    // count total nodes on edges
    int f,totalNode = 0;
    for(f=0;f<Nedge;f++){
        int f1 = edge->varfM[f];
        totalNode += dg_cell_Nfp(cell)[f1];
    }
    // map node
    int *varpM = vector_int_create(totalNode);
    int *varpP = vector_int_create(totalNode);
    int *varfpM = vector_int_create(totalNode);
    int *varfpP = vector_int_create(totalNode);

    int *Nfpstart = (int *)calloc(Nfaces, sizeof(int));
    Nfpstart[0] = 0;
    for(f=1;f<Nfaces;f++){
        Nfpstart[f] = Nfpstart[f-1] + dg_cell_Nfp(cell)[f-1];
    }

    int n1,n2,sk=0;
    for(f=0;f<Nedge;f++){
        int k1 = edge->varkM[f];
        int k2 = edge->varkP[f];
        int f1 = edge->varfM[f];
        int f2 = edge->varfP[f];
        int ftype = edge->ftype[f];

        const int Nfp = dg_cell_Nfp(cell)[f1];
        for(n1=0;n1<Nfp;n1++){
            const int idM = k1*Np + Fmask[f1][n1];
            varpM[sk] = idM;
            varfpM[sk] = Nfpstart[f1] + n1;
            double xM = x[0][idM];
            double yM = y[0][idM];
            double zM = z[0][idM];

            if(ftype != FACE_PARALL){
                for(n2=0;n2<Nfp;n2++){
                    const int idP = k2*Np + Fmask[f2][n2];
                    double xP = x[0][idP];
                    double yP = y[0][idP];
                    double zP = z[0][idP];

                    double d12 = (xM-xP)*(xM-xP)+(yM-yP)*(yM-yP)+(zM-zP)*(zM-zP);
                    if(d12 < EPS){
                        varpP[sk] = idP;
                        varfpP[sk] = Nfpstart[f2] + n2;
                        break;
                    }
                }
            }else{ // parallel face
                for(n2=0;n2<TotalParNode;n2++){
                    if( varpM[sk] == parnode[n2] ){
                        varpP[sk] = n2;
                        varfpP[sk] = 0;
                        break;
                    }
                }
            }
            sk++;
        }
    }
    // assignment
    edge->Nnode = totalNode;
    edge->varfpM = varfpM;
    edge->varfpP = varfpP;
    edge->varpM = varpM;
    edge->varpP = varpP;

    return;
}