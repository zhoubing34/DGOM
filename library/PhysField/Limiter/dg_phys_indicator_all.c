//
// Created by li12242 on 17/4/20.
//

#include "dg_phys_indicator_all.h"

void dg_phys_indicator_all(dg_phys_info *phys, int *tind){
    dg_grid *grid = dg_phys_info_grid(phys);
    const int K = dg_grid_K(grid);
    const int Nfield = phys->Nfield;
    const int Nfaces = dg_cell_Nfaces(dg_phys_info_cell(phys));

    int k,fld,sk=0;
    for(k=0;k<K;k++){
        for(fld=0;fld<Nfield;fld++){
            tind[sk++] = 1;
        }
    }
    return;
}