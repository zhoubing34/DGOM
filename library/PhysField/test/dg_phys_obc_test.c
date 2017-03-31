//
// Created by li12242 on 17/2/23.
//

#include "dg_phys_obc_test.h"

int dg_phys_obc_test(dg_phys *phys, int verbose){
    dg_cell *cell = dg_phys_cell(phys);
    const int procid = dg_mesh_procid(dg_phys_mesh(phys));
    const int nprocs = dg_mesh_nprocs(dg_phys_mesh(phys));
    const int Nfield = dg_phys_Nfield(phys);
    const int K = dg_grid_K(dg_phys_grid(phys));
    const int Np = dg_cell_Np(cell);
    const int Nfaces = dg_cell_Nfaces(cell);
    int fail=0;

    /* get the cell and face */
    int **EToBS = dg_grid_EToBS(dg_phys_grid(phys));
    int k,n,f,fld,NOBnodes=0;
    for(k=0;k<K;k++){
        for(f=0;f<Nfaces;f++){
            if(EToBS[k][f] >= FACE_OPENBS){
                int Nfp = dg_cell_Nfp(cell)[f];
                NOBnodes += Nfp;
            }
        }
    }
    int *OBnodes = vector_int_create(NOBnodes);
    int sk=0;
    for(k=0;k<K;k++){
        for(f=0;f<Nfaces;f++){
            if(EToBS[k][f] >= FACE_OPENBS){
                int Nfp = dg_cell_Nfp(cell)[f];
                for(n=0;n<Nfp;n++){
                    OBnodes[sk++] = k*Np + dg_cell_Fmask(cell)[f][n];
                }
            }
        }
    }

    double ftime = 44714, dt = 600, time;
    /* print the open boundary value in files */
    FILE *fp = create_log(__FUNCTION__, procid, nprocs);
    /* print coordinates */
    fprintf(fp, "x = [");
    for(n=0;n<NOBnodes;n++){
        fprintf(fp, "%f, ", dg_phys_region(phys)->x[0][OBnodes[n]]);
    }
    fprintf(fp, "]\n\ny = [");
    for(n=0;n<NOBnodes;n++){
        fprintf(fp, "%f, ", dg_phys_region(phys)->y[0][OBnodes[n]]);
    }
    fprintf(fp, "]\n\n");

    /* print values */
    for(time=-0.1;time<ftime;time+=dt){
        phys->obc_update(phys, time);
        fprintf(fp, "time = %f\nf_ext = [", time);
        for(n=0;n<NOBnodes;n++){
            for(fld=0;fld<Nfield;fld++){
                fprintf(fp, "%f, ", dg_phys_f_extQ(phys)[OBnodes[n]*Nfield + fld]);
            }
        }
        fprintf(fp, "]\n\n");
    }
    fclose(fp);

    if(!procid){
        if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    }
    return fail;
}