//
// Created by li12242 on 17/1/5.
//

#include "dg_phys_test.h"
#include "../dg_phys_strong_LDG_opt.h"

static int wall_func(dg_real nx, dg_real ny, dg_real *varM, dg_real *varP){
    varP[0] = varM[0];
    varP[1] = varM[1];
    return 0;
}

static int obc_func(dg_real nx, dg_real ny, dg_real *f_M, dg_real *f_ext, int obc_ind, dg_real *f_P){
    f_P[0] = f_M[0];
    f_P[1] = f_M[1];
    return 0;
}

static int vis_func(dg_real *f_Q, dg_real *vis_Q){
    vis_Q[0] = f_Q[0];
    vis_Q[1] = f_Q[1];
    return 0;
}

int dg_phys_strong_LDG_opt2d_test(dg_phys *phys, int verbose){
    int fail = 0;
    dg_mesh *mesh = dg_phys_mesh(phys);
    dg_region *region = dg_phys_region(phys);
    dg_cell *cell = dg_phys_cell(phys);

    const int Nfield = dg_phys_Nfield(phys);
    const int K = dg_grid_K(dg_phys_grid(phys));
    const int Np = dg_cell_Np(cell);
    const int procid = dg_region_procid(region);
    const int nprocs = dg_region_nprocs(region);

    int k,n;
    const dg_real miu = 0.5;
    const dg_real sqrt_miu = dg_sqrt(miu);
    dg_real px_ext[Np*Nfield*K], py_ext[Np*Nfield*K], rhs_ext[Np*Nfield*K];

    dg_real *f_rhsQ = dg_phys_f_rhsQ(phys);
    dg_real *f_Q = dg_phys_f_Q(phys);

    dg_phys_LDG *ldg = phys->ldg;
    dg_real *miu_x = dg_phys_ldg_sqrt_miux(ldg);
    dg_real *miu_y = dg_phys_ldg_sqrt_miuy(ldg);
    // initialize
    for(k=0;k<K;k++){
        for(n=0;n<Np;n++){
            const int sk = k*Np*Nfield + n*Nfield;

            dg_real xt = region->x[k][n];
            dg_real yt = region->y[k][n];

            f_Q[sk] = xt*xt;
            f_rhsQ[sk] = 0;
            miu_x[sk] = sqrt_miu; miu_y[sk] = sqrt_miu;
            px_ext[sk] = 2.*xt*sqrt_miu; py_ext[sk] = 0.0*sqrt_miu;
            rhs_ext[sk] = 2.0*miu;

            const int st = sk+1;
            f_rhsQ[st] = 0;
            f_Q[st] = xt*yt+yt*yt;
            miu_x[st] = sqrt_miu; miu_y[st] = sqrt_miu;
            px_ext[st] = yt*sqrt_miu; py_ext[st] = (xt+2.*yt)*sqrt_miu;
            rhs_ext[st] = 2.0*miu;
        }
    }

    MPI_Request mpi_send_requests[nprocs];
    MPI_Request mpi_recv_requests[nprocs];
    int Nmess=0;
    /* fetch nodal value through all procss */
    Nmess = phys->fetch_node_buffer(phys, mpi_send_requests, mpi_recv_requests);
    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_recv_requests, instatus);
    MPI_Waitall(Nmess, mpi_send_requests, instatus);

    dg_phys_LDG_solve_vis_opt2d(phys, vis_func, wall_func, wall_func, obc_func);
    //dg_phys_LDG_solve_vis_opt2d(phys, vis_func, wall_func, wall_func, obc_func);

    if(!procid){
        fail = vector_double_test(__FUNCTION__, dg_phys_ldg_px(ldg), px_ext, Np * Nfield * K);
        fail = vector_double_test(__FUNCTION__, dg_phys_ldg_py(ldg), py_ext, Np * Nfield * K);
        fail = vector_double_test(__FUNCTION__, dg_phys_f_rhsQ(phys), rhs_ext, Np * Nfield * K);
    }

    if(verbose){
        FILE *fp = create_log(__FUNCTION__, mesh->procid, mesh->nprocs);
        fprintf(fp, "K = %d\n", K);
        fprintf(fp, "Nfield = %d\n", Nfield);
        fprintf(fp, "Np = %d\n", Np);
        print_double_vector2file(fp, "f_Q", f_Q, K*Np*Nfield);
        print_double_vector2file(fp, "px_ext", px_ext, Nfield * Np * K);
        print_double_vector2file(fp, "px_Q", dg_phys_ldg_px(ldg), Nfield * Np * K);
        print_double_vector2file(fp, "py_ext", py_ext, Nfield * Np * K);
        print_double_vector2file(fp, "py_Q", dg_phys_ldg_px(ldg), Nfield * Np * K);
        print_double_vector2file(fp, "rhsQ_ext", rhs_ext, Nfield * Np * k);
        print_double_vector2file(fp, "f_rhsQ", dg_phys_ldg_px(ldg), Nfield * Np * K);
        fclose(fp);
    }

    if(!procid) {
        if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    }
    return fail;
}