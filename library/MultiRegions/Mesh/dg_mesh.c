//
// Created by li12242 on 12/19/16.
//

#include "dg_mesh.h"
#include "dg_mesh_connect.h"
#include "dg_mesh_fetch_buffer.h"

static void dg_mesh_free2d(dg_mesh *mesh);

typedef struct dg_mesh_creator{
    void (*init_cell_fetch_buffer)(dg_mesh *mesh);
    void (*init_node_fetch_buffer)(dg_mesh *mesh);
    int (*fetch_node_buffer)(dg_mesh *mesh, int Nfield, dg_real *f_Q, dg_real *f_recvQ,
                             MPI_Request *send_requests,
                             MPI_Request *recv_requests);
    int (*fetch_cell_buffer)(dg_mesh *mesh, int Nfield, dg_real *f_Q, dg_real *f_recvQ,
                             MPI_Request *send_requests,
                             MPI_Request *recv_requests);
    void (*free_func)(dg_mesh *mesh);
}dg_mesh_creator;

const static dg_mesh_creator mesh2d_creator = {
        dg_mesh_init_cell_fetch_buffer,
        dg_mesh_init_node_fetch_buffer,
        dg_mesh_fetch_node_buffer,
        dg_mesh_fetch_cell_buffer,
        dg_mesh_free2d,
};

dg_mesh* dg_mesh_create(dg_region *region){

    dg_mesh *mesh = (dg_mesh *)calloc(1, sizeof(dg_mesh));
    /* basic information */
    mesh->nprocs = region->nprocs;
    mesh->procid = region->procid;
    mesh->region = region;
    mesh->grid = region->grid;
    mesh->cell = region->cell;

    dg_cell_type type = dg_cell_celltype(region->cell);
    const dg_mesh_creator *creator;
    switch (type){
        case TRIANGLE:
            creator = &mesh2d_creator; break;
        case QUADRIL:
            creator = &mesh2d_creator; break;
        default:
            fprintf(stderr, "%s (%d):\nUnknown cell type %d\n", __FUNCTION__, __LINE__, type);
            exit(-1);
    }
    creator->init_cell_fetch_buffer(mesh);
    creator->init_node_fetch_buffer(mesh);
    mesh->fetch_node_buffer = creator->fetch_node_buffer;
    mesh->fetch_cell_buffer = creator->fetch_cell_buffer;
    mesh->free_func = creator->free_func;
    return mesh;
}

static void dg_mesh_free2d(dg_mesh *mesh){
    /* cell connection */
    vector_int_free(mesh->parfaceNum);
    vector_int_free(mesh->parcell);
    vector_int_free(mesh->parface);
    /* node connection */
    vector_int_free(mesh->parnodeNum);
    vector_int_free(mesh->parnode);
    free(mesh);
    return;
}

void dg_mesh_free(dg_mesh *mesh){
    mesh->free_func(mesh);
    return;
}