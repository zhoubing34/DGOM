//
// Created by li12242 on 17/3/10.
//

#include "dg_edge_test.h"

int dg_edge_facemap_test(dg_edge *edge, int verbose){
    int fail = 0;
    int procid,nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    if(verbose){
        FILE *fp = create_log(__FUNCTION__, nprocs, procid);
        fprintf(fp, "Nedge = %d\n", edge->Nedge);
        print_int_vector2file(fp, "varkM", edge->varkM, edge->Nedge);
        print_int_vector2file(fp, "varkP", edge->varkP, edge->Nedge);
        print_int_vector2file(fp, "varfM", edge->varfM, edge->Nedge);
        print_int_vector2file(fp, "varfP", edge->varfP, edge->Nedge);
        print_int_vector2file(fp, "ftype", edge->ftype, edge->Nedge);
        fclose(fp);
    }
    if( (!procid) & (!fail) ) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_edge_nodemap_test(dg_edge *edge, int verbose){
    int fail = 0;
    int procid,nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    double **x = dg_region_x(dg_edge_region(edge));
    double **y = dg_region_y(dg_edge_region(edge));
    const int Nnode = dg_edge_Nnode(edge);
    const int Nedge = dg_edge_Nedge(edge);
    double xM[Nnode], xP[Nnode];
    double yM[Nnode], yP[Nnode];
    int n,f,sk=0;
    for(f=0;f<Nedge;f++){
        int f1 = edge->varfM[f];
        int ftype = edge->ftype[f];
        const int Nfp = dg_cell_Nfp(dg_edge_cell(edge))[f1];
        for(n=0;n<Nfp;n++){
            if(ftype != FACE_PARALL){
                xM[sk] = x[0][edge->varpM[sk]];
                xP[sk] = x[0][edge->varpP[sk]];
                yM[sk] = y[0][edge->varpM[sk]];
                yP[sk] = y[0][edge->varpP[sk]];
            }else{
                xM[sk] = x[0][edge->varpM[sk]];
                xP[sk] = x[0][edge->varpM[sk]];
                yM[sk] = y[0][edge->varpM[sk]];
                yP[sk] = y[0][edge->varpM[sk]];
            }
            sk++;
        }
    }

    fail = vector_double_test(__FUNCTION__, xP, xM, Nnode);
    fail = vector_double_test(__FUNCTION__, yP, yM, Nnode);

    if(verbose){
        FILE *fp = create_log(__FUNCTION__, nprocs, procid);
        fprintf(fp, "Nnode = %d\n", Nnode);
        print_int_vector2file(fp, "varpM", edge->varpM, Nnode);
        print_int_vector2file(fp, "varpP", edge->varpP, Nnode);
        fclose(fp);
    }
    if( (!procid) & (!fail) ) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_edge_surfinfo_test(dg_edge *edge, int verbose){
    int fail = 0;
    int procid,nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    if(verbose){
        FILE *fp = create_log(__FUNCTION__, nprocs, procid);
        fprintf(fp, "Nnode = %d\n", edge->Nnode);
        print_double_vector2file(fp, "nx", edge->nx, edge->Nnode);
        print_double_vector2file(fp, "ny", edge->ny, edge->Nnode);
        print_double_vector2file(fp, "fsc", edge->fsc, edge->Nnode);
        fclose(fp);
    }
    if( (!procid) & (!fail) ) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}