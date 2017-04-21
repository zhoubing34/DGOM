//
// Created by li12242 on 17/3/10.
//

#include "dg_edge.h"
#include "dg_edge_face_map.h"
#include "dg_edge_node_map.h"
#include "dg_edge_surfinfo.h"

static void dg_edge_set_info2d(dg_edge *edge);

/**
 * @brief creator for generating dg_edge structure.
 */
typedef struct dg_edge_creator{
    void (*map_face)(dg_edge *edge);
    void (*map_node)(dg_edge *edge);
    void (*get_surf_info)(dg_edge *edge);
    void (*set_info)(dg_edge *edge);
}dg_edge_creator;

static const dg_edge_creator edge_creator2d = {
        dg_edge_face_map2d,
        dg_edge_node_map,
        dg_edge_surfinfo2d,
        dg_edge_set_info2d,
};

dg_edge* dg_edge_create(dg_mesh *mesh){
    dg_edge *edge = (dg_edge *) calloc(1, sizeof(dg_edge));
    /* basic info */
    edge->mesh = mesh;

    const dg_edge_creator *creator;
    dg_cell_type type = dg_cell_celltype(dg_edge_cell(edge));
    switch (type){
        case TRIANGLE: creator = &edge_creator2d; break;
        case QUADRIL: creator = &edge_creator2d; break;
        default:
            fprintf(stderr, "%s (%d): Unknown cell type %d\n",
                    __FUNCTION__, __LINE__, type);
            exit(-1);
    }
    creator->map_face(edge);
    creator->map_node(edge);
    creator->get_surf_info(edge);
    creator->set_info(edge);
    return edge;
}

static void dg_edge_set_info2d(dg_edge *edge){
    const int Nedge = dg_edge_Nedge(edge);
    const int Nnode = dg_edge_Nnode(edge);
    int *faceinfo = (int*) calloc(Nedge*5, sizeof(int));
    dg_real *nodeinfo = (dg_real*) calloc(Nnode*7, sizeof(dg_real));

    int f,n,sk=0;
    for(f=0;f<Nedge;f++){
        faceinfo[sk++] = edge->varkM[f];
        faceinfo[sk++] = edge->varkP[f];
        faceinfo[sk++] = edge->varfM[f];
        faceinfo[sk++] = edge->varfP[f];
        faceinfo[sk++] = edge->ftype[f];
    }
    sk = 0;
    for(n=0;n<Nnode;n++){
        nodeinfo[sk++] = edge->varpM[n];
        nodeinfo[sk++] = edge->varpP[n];
        nodeinfo[sk++] = edge->varfpM[n];
        nodeinfo[sk++] = edge->varfpP[n];
        nodeinfo[sk++] = edge->nx[n];
        nodeinfo[sk++] = edge->ny[n];
        nodeinfo[sk++] = edge->fsc[n];
    }

    edge->surfinfo = faceinfo;
    edge->nodeinfo = nodeinfo;
    return;
}

void dg_edge_free(dg_edge *edge){
    vector_int_free(edge->varkM);
    vector_int_free(edge->varkP);
    vector_int_free(edge->varfM);
    vector_int_free(edge->varfP);
    vector_int_free(edge->ftype);

    vector_int_free(edge->varpM);
    vector_int_free(edge->varpP);

    vector_real_free(edge->fsc);
    vector_real_free(edge->nx);
    vector_real_free(edge->ny);
    vector_real_free(edge->nz);
    free(edge);
    return;
}
