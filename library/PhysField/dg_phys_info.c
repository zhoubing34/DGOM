//
// Created by li12242 on 17/3/22.
//

#include "dg_phys_info.h"

static int dg_phys_fetch_node_buffer(dg_phys_info *phys_info, MPI_Request *mpi_send_requests,
                                     MPI_Request *mpi_recv_requests);
static int dg_phys_fetch_cell_buffer(dg_phys_info *phys_info, MPI_Request *mpi_send_requests,
                                     MPI_Request *mpi_recv_requests);
static void dg_phys_info_cell_mean(dg_phys_info *phys);

/**
 * @brief create a new pointer to dg_phys_info structure.
 * @param Nfields number of physical field;
 * @param edge pointer to a dg_edge structure;
 * @return
 * phys_info pointer to a dg_phys_info structure.
 */
dg_phys_info* dg_phys_info_create(int Nfields, dg_area *area){
    dg_phys_info *phys_info = (dg_phys_info*)calloc(1, sizeof(dg_phys_info));

    phys_info->area = area;
    phys_info->Nfield = Nfields; ///< number of fields

    const int K = dg_grid_K(dg_area_grid(area));
    const int Np = dg_cell_Np(dg_area_cell(area));
    const int NfetchFace = dg_mesh_NfetchFace(dg_area_mesh(area));
    const int NfetchNode = dg_mesh_NfetchNode(dg_area_mesh(area));

    /* nodal array */
    phys_info->f_Q    = (dg_real *) calloc((size_t) K*Np*Nfields, sizeof(dg_real));
    phys_info->f_recvQ = (dg_real *) calloc((size_t) NfetchNode * Nfields, sizeof(dg_real));

    /* volume array */
    phys_info->c_Q = (dg_real *) calloc((size_t) K*Nfields, sizeof(dg_real));
    phys_info->c_recvQ = (dg_real *) calloc((size_t) NfetchFace * Nfields, sizeof(dg_real));

    /* time discretize */
    phys_info->f_rhsQ = (dg_real *) calloc((size_t) K*Np*Nfields, sizeof(dg_real));
    phys_info->f_resQ = (dg_real *) calloc((size_t) K*Np*Nfields, sizeof(dg_real));

    /* methods */
    phys_info->fetch_node_buffer = dg_phys_fetch_node_buffer;
    phys_info->fetch_cell_buffer = dg_phys_fetch_cell_buffer;
    phys_info->cell_mean = dg_phys_info_cell_mean;

    return phys_info;
}

void dg_phys_info_free(dg_phys_info *phys_info){
    free(phys_info->f_Q);
    free(phys_info->f_recvQ);

    free(phys_info->c_Q);
    free(phys_info->c_recvQ);

    free(phys_info->f_rhsQ);
    free(phys_info->f_resQ);

    free(phys_info);
    return;
}

/**
 * @brief Send/rece nodal value node information through buffers
 *
 * @details
 * The internal boundary values of `f_Q` in local process is arranged into `f_outQ` buffer
 * to send to other processes. While `mpi_send_requests` get the request of each
 * MPI_Isend function. The incoming boundary values is stored into the `f_inQ` buffer
 * and the request are stored in `mpi_recv_requests`.
 * `Nmess` gets the number of send and receive requests.
 *
 * Usages:
 *
 *     MPI_Request mpi_out_requests[nprocs];
 *     MPI_Request mpi_in_requests[nprocs];
 *     int Nmess;
 *     Nmess = phys_fetchNodeBuffer2d(phys, mpi_out_requests, mpi_in_requests);
 *
 *     MPI_Status instatus[nprocs];
 *     MPI_Waitall(Nmess, mpi_in_requests, instatus);
 *
 * @param[in] phys_info pointer to a dg_phys_info structure;
 * @param[in,out] mpi_send_requests MPI_Request pointer;
 * @param[in,out] mpi_recv_requests MPI_Request pointer;
 * @return
 * Nmess number of messages in `mpi_send_requests` and `mpi_send_requests`
 *
 */
static int dg_phys_fetch_node_buffer(dg_phys_info *phys_info,
                                     MPI_Request *mpi_send_requests,
                                     MPI_Request *mpi_recv_requests) {
    dg_mesh *mesh = dg_phys_info_mesh(phys_info);
    const int Nfield = phys_info->Nfield;
    int Nmess;
    Nmess = mesh->fetch_node_buffer(mesh, Nfield, phys_info->f_Q, phys_info->f_recvQ,
                                    mpi_send_requests, mpi_recv_requests);
    return Nmess;
}

/**
 * @brief Send/rece elemental value `c_Q` through buffers
 *
 * @details
 * The internal boundary elemental values of `c_Q` in local process is arranged into
 * `c_outQ` buffer to send to other processes. While `mpi_send_requests` get
 * the request of each `MPI_Isend` function. The incoming elemental values is stored
 * into the `c_inQ` buffer and the request are stored in `mpi_recv_requests`.
 * `Nmess` gets the number of send and receive requests.
 *
 * Usages:
 *
 *     MPI_Request mpi_out_requests[nprocs];
 *     MPI_Request mpi_in_requests[nprocs];
 *     int Nmess;
 *     phys_fetchCellBuffer(phys, mpi_out_requests, mpi_in_requests, &Nmess);
 *
 *     MPI_Status instatus[nprocs];
 *     MPI_Waitall(Nmess, mpi_in_requests, instatus);
 *
 * @param[in] phys_info pointer to a dg_phys_info structure;
 * @param[in,out] mpi_send_requests MPI_Request pointer;
 * @param[in,out] mpi_recv_requests MPI_Request pointer;
 * @return Nmess number of messages in `mpi_send_requests` and `mpi_send_requests`.
 *
 */
static int dg_phys_fetch_cell_buffer(dg_phys_info *phys_info,
                                     MPI_Request *mpi_send_requests,
                                     MPI_Request *mpi_recv_requests){

    dg_mesh *mesh = dg_phys_info_mesh(phys_info);
    const int Nfield = phys_info->Nfield;
    int Nmess;
    Nmess = mesh->fetch_cell_buffer(mesh, Nfield, phys_info->c_Q, phys_info->c_recvQ,
                                    mpi_send_requests, mpi_recv_requests);
    return Nmess;
}


/**
 * @brief calculate the integral averaged cell value of each physical field.
 * @param[in,out] phys pointer to dg_phys structure;
 */
static void dg_phys_info_cell_mean(dg_phys_info *phys){

    dg_region *region = dg_phys_info_region(phys);
    const int Nfield = phys->Nfield;
    const int K = dg_grid_K(dg_phys_info_grid(phys));
    const int Np = dg_cell_Np(dg_phys_info_cell(phys));

    register int k,fld;

    dg_real *f_Q = phys->f_Q;
    dg_real *c_Q = phys->c_Q;

    for(k=0;k<K;k++){
        region->vol_integral(region, Nfield, k, f_Q+k*Np*Nfield, c_Q+k*Nfield);
        const double Area = 1.0/region->size[k];
        for(fld=0;fld<Nfield;fld++){ // divided the area
            c_Q[k*Nfield + fld] *= Area;
        }
    }

    return;
}
