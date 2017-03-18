//
// Created by li12242 on 17/3/9.
//
#include "dg_grid.h"
#include "dg_grid_BS.h"

/**
 * @brief Initialize the EToBS field.
 * @param grid dg_grid structure;
 */
void dg_grid_init_BS(dg_grid *grid){
    const int procid = dg_grid_procid(grid);
    const int K = dg_grid_K(grid);
    const int Nfaces = dg_cell_Nfaces(grid->cell);

    int **EToBS = matrix_int_create(K, Nfaces);
    int **EToP = grid->EToP;
    int k,f;
    for(k=0;k<K;k++){
        for(f=0;f<Nfaces;f++){
            /* the inner boundary surface */
            if(EToP[k][f] != procid){
                EToBS[k][f] = FACE_PARALL;
            }else{
                EToBS[k][f] = FACE_INNER;
            }
        }
    }
    grid->EToBS = EToBS;
    return;
}

void dg_grid_add_BS3d(dg_grid *grid, int Nsurface, int **SFToV){
    return;
}
/**
 * @brief Add boundary surface type into EToBS.
 * @details
 *
 * @param grid pointer to dg_grid structure;
 * @param Nsurface number of surface in SFToV;
 * @param SFToV surface type to vertex, size [Nsurf x 3];
 */
void dg_grid_add_BS2d(dg_grid *grid, int Nsurface, int **SFToV){

    const int K = dg_grid_K(grid);
    const int Nvert = dg_grid_Nv(grid);
    const int Nfaces = dg_cell_Nfaces(grid->cell);
    const int Nv = dg_cell_Nv(grid->cell);

    int **EToV = grid->EToV;
    int **EToBS = grid->EToBS;
    /* check and sort the vertex sequence in SFToV */
    int k,f;
    int face_indicator[Nsurface];
    for(f=0;f<Nsurface;f++){
        int v1 = SFToV[f][0];
        int v2 = SFToV[f][1];
        face_indicator[f] = min(v1, v2)*Nvert + max(v1, v2);
        int facetype = SFToV[f][2];
        if( facetype == FACE_INNER | facetype == FACE_PARALL ){
            printf("%s (%d)\nWarning: surface type in SFToV[%d][2] = %d\n",
                   __FUNCTION__, __LINE__, f, SFToV[f][2]);
        }
    }

    // pair the face in EToBS
    int f1;
    for(k=0;k<K;k++){
        for(f=0;f<Nfaces;f++){
            /* get the usr-defined boundary */
            int n1 = f;
            int n2 = (f+1)%Nv;
            int t1 = EToV[k][n1];
            int t2 = EToV[k][n2];
            int t_temp = min(t1, t2)*Nvert + max(t1, t2);
            for(f1=0;f1<Nsurface;f1++){
                if(t_temp == face_indicator[f1]){ // compare the face2d
                    EToBS[k][f] = SFToV[f1][2];
                    break; // jump out loop of SFToV
                }
            }
        }
    }

    return;
}