//
// Created by li12242 on 17/1/1.
//

#include "dg_phys_strong_LDG_opt.h"

#define DEBUG 0
//#define DELTA_1 0
//#define DELTA_2 0
//#define DELTA_3 0

#if DEBUG
#include "LibUtilities/UTest.h"
#endif

static void ldg_auxi_vol_opt2d(dg_phys *phys, dg_phys_LDG *ldg, Vis_Fun vis_func);
static void ldg_auxi_surf_opt2d(dg_phys *phys, dg_phys_LDG *ldg, Vis_Fun vis_func,
                                Wall_Condition slipwall_func, Wall_Condition nonslipwall_func, OBC_Fun obc_fun);
static void ldg_phys_vol_opt2d(dg_phys *phys, dg_phys_LDG *ldg);
static void ldg_auxi_times_vis2d(dg_phys *phys, dg_phys_LDG *ldg);
static void ldg_phys_sruf_opt2d(dg_phys *phys, dg_phys_LDG *ldg);

void dg_phys_LDG_solve_vis_opt2d(dg_phys *phys, Vis_Fun vis_func,
                                 Wall_Condition slipwall_func,
                                 Wall_Condition non_slipwall_func,
                                 OBC_Fun obc_fun){

    dg_mesh *mesh = dg_phys_mesh(phys);
    dg_phys_LDG *ldg = phys->ldg;
    const int nprocs = dg_phys_nprocs(phys);

    /* 0. send/recv the viscosity value */
    const int Nfield = dg_phys_Nfield(phys);
    int Nmess;
    MPI_Request mpi_x_send_requests[nprocs], mpi_x_recv_requests[nprocs];
    MPI_Request mpi_y_send_requests[nprocs], mpi_y_recv_requests[nprocs];
    /* 1. calculate the auxiliary variables - volume integral */
    ldg_auxi_vol_opt2d(phys, ldg, vis_func);

    /* 2. calculate the auxiliary variables - surface integral */
    ldg_auxi_surf_opt2d(phys, ldg, vis_func, slipwall_func, non_slipwall_func, obc_fun);

    /* 3. multiply sqrt root of viscosity to auxiliary variables */
    ldg_auxi_times_vis2d(phys, ldg);

    /* 4. send/recv auxiliary variables */
    Nmess = mesh->fetch_node_buffer(mesh, Nfield, dg_phys_ldg_px(ldg), ldg->px_recv,
                                    mpi_x_send_requests, mpi_x_recv_requests);
    Nmess = mesh->fetch_node_buffer(mesh, Nfield, dg_phys_ldg_py(ldg), ldg->py_recv,
                                    mpi_y_send_requests, mpi_y_recv_requests);

    /* 5. calculate the volume term */
    ldg_phys_vol_opt2d(phys, ldg);

    /* wait for message passing */
    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_x_recv_requests, instatus);
    MPI_Waitall(Nmess, mpi_y_recv_requests, instatus);

    /* 6. calculate the surface term */
    ldg_phys_sruf_opt2d(phys, ldg);

    /* wait for message passing */
    MPI_Waitall(Nmess, mpi_x_send_requests, instatus);
    MPI_Waitall(Nmess, mpi_y_send_requests, instatus);
    return;
}

static void ldg_auxi_times_vis2d(dg_phys *phys, dg_phys_LDG *ldg){
    dg_cell *cell = dg_phys_cell(phys);
    dg_grid *grid = dg_phys_grid(phys);

    const int K = dg_grid_K(grid);
    const int Np = dg_cell_Np(cell);
    const int Nfield = dg_phys_Nfield(phys);

    int k,n,fld,sk=0;
    for(k=0;k<K;k++){
        for(n=0;n<Np;n++){
            for(fld=0;fld<Nfield;fld++){
                dg_real miu_x = dg_phys_ldg_sqrt_miux(ldg)[sk];
                dg_real miu_y = dg_phys_ldg_sqrt_miuy(ldg)[sk];
                dg_phys_ldg_px(ldg)[sk] *= miu_x;
                dg_phys_ldg_py(ldg)[sk] *= miu_y;
                sk++;
            }
        }
    }
    return;
}

static void ldg_auxi_vol_opt2d(dg_phys *phys, dg_phys_LDG *ldg, Vis_Fun vis_func){

    dg_cell *cell = dg_phys_cell(phys);
    dg_grid *grid = dg_phys_grid(phys);
    dg_region *region = dg_phys_region(phys);

    const int K = dg_grid_K(grid);
    const int Np = dg_cell_Np(cell);
    const int Nfield = dg_phys_Nfield(phys);

    dg_real *f_Q = dg_phys_f_Q(phys); //phys->f_Q;
    dg_real *f_Dr = dg_cell_f_Dr(cell); //phys->cell->f_Dr;
    dg_real *f_Ds = dg_cell_f_Ds(cell); //phys->cell->f_Ds;
    dg_real *px = dg_phys_ldg_px(ldg);
    dg_real *py = dg_phys_ldg_py(ldg);

    register unsigned int k,n,m,fld;
    dg_real **drdx_p = region->drdx;
    dg_real **drdy_p = region->drdy;
    dg_real **dsdx_p = region->dsdx;
    dg_real **dsdy_p = region->dsdy;

    dg_real *vis_term = vector_real_create(Np*Nfield);
    dg_real *px_rhs = vector_real_create(Nfield);
    dg_real *py_rhs = vector_real_create(Nfield);

    int rhsid = 0;
    for(k=0;k<K;k++){
        dg_real *var = f_Q + k*Np*Nfield; // variable in k-th element
        // calculate viscosity term
        for(n=0;n<Np;n++){
            vis_func(var+n*Nfield, vis_term+n*Nfield);
        }

        for(n=0;n<Np;n++){
            const dg_real *ptDr = f_Dr+n*Np; // n-th row of Dr
            const dg_real *ptDs = f_Ds+n*Np; // n-th row of Ds

            const dg_real drdx = drdx_p[k][n]; // volume geometry for n-th point
            const dg_real drdy = drdy_p[k][n]; // volume geometry for n-th point
            const dg_real dsdx = dsdx_p[k][n]; // volume geometry for n-th point
            const dg_real dsdy = dsdy_p[k][n]; // volume geometry for n-th point

            // initialize rhs
            for(m=0;m<Nfield;m++){
                px_rhs[m] = 0;
                py_rhs[m] = 0;
            }

            for(m=0;m<Np;++m){
                const dg_real dr = ptDr[m]; // m-th column for Dr
                const dg_real ds = ptDs[m]; // m-th column for Ds
                const dg_real dx = drdx*dr+dsdx*ds;
                const dg_real dy = drdy*dr+dsdy*ds;

                const dg_real *vis = vis_term + m*Nfield;

                for(fld=0;fld<Nfield;fld++){
                    px_rhs[fld] += dx*vis[fld];
                    py_rhs[fld] += dy*vis[fld];
                }
            }

            for(fld=0;fld<Nfield;fld++){
                px[ rhsid   ] = px_rhs[fld];
                py[ rhsid++ ] = py_rhs[fld];
            }

        }

    }

    vector_real_free(vis_term);
    vector_real_free(px_rhs);
    vector_real_free(py_rhs);

    return;
}

static void ldg_auxi_surf_opt2d(dg_phys *phys, dg_phys_LDG *ldg, Vis_Fun vis_func,
                                Wall_Condition slipwall_func,
                                Wall_Condition non_slipwall_func,
                                OBC_Fun obc_fun){

    dg_edge *edge = dg_phys_edge(phys);
    dg_cell *cell = dg_phys_cell(phys);
    const int Nfield = dg_phys_Nfield(phys);
    const int Nedge = dg_edge_Nedge(edge);
    const int Nfptotal = dg_cell_Nfptotal(cell);
    const int Np = dg_cell_Np(cell);

    dg_real *f_Q = dg_phys_f_Q(phys);
    dg_real *f_inQ = dg_phys_f_recvQ(phys); // phys->f_recvQ;
    dg_real *f_ext = dg_phys_f_extQ(phys); // phys->f_extQ;
    dg_real *f_LIFT  = dg_cell_f_LIFT(cell); //phys->cell->f_LIFT;

    dg_real *px = dg_phys_ldg_px(ldg);
    dg_real *py = dg_phys_ldg_py(ldg);

    register int f,m,n,fld,surfid=0,nodeid=0;

    for(f=0;f<Nedge;f++){
        const int k1 = edge->surfinfo[surfid++];
        const int k2 = edge->surfinfo[surfid++];
        const int f1 = edge->surfinfo[surfid++];
        const int f2 = edge->surfinfo[surfid++];
        const int ftype = edge->surfinfo[surfid++];
        const int Nfp = dg_cell_Nfp(cell)[f1];

        dg_real pxFlux_M[Nfp*Nfield], pxFlux_P[Nfp*Nfield];
        dg_real pyFlux_M[Nfp*Nfield], pyFlux_P[Nfp*Nfield];
        int fp_M[Nfp], fp_P[Nfp];
        for(m=0;m<Nfp;m++){
            const int idM = (int)edge->nodeinfo[nodeid++];
            const int idP = (int)edge->nodeinfo[nodeid++];
            fp_M[m] = (int)edge->nodeinfo[nodeid++];
            fp_P[m] = (int)edge->nodeinfo[nodeid++];
            const dg_real nx = edge->nodeinfo[nodeid++];
            const dg_real ny = edge->nodeinfo[nodeid++];
            const dg_real fsc = edge->nodeinfo[nodeid++];

            dg_real f_M[Nfield], f_P[Nfield];
            // local face2d values
            for(fld=0;fld<Nfield;fld++) { f_M[fld] = f_Q[idM*Nfield+fld]; }
            // adjacent nodal values
            switch (ftype){
                case FACE_INNER:
                    for(fld=0;fld<Nfield;fld++){ f_P[fld] = f_Q[idP*Nfield+fld]; }
                    break;
                case FACE_PARALL:
                    for(fld=0;fld<Nfield;fld++){ f_P[fld] = f_inQ[idP*Nfield+fld]; }
                    break;
                case FACE_SLIPWALL:
                    slipwall_func(nx, ny, f_M, f_P);
                    break;
                case FACE_NSLIPWALL:
                    non_slipwall_func(nx, ny, f_M, f_P);
                    break;
                default: /// open boundary condition
                    obc_fun(nx, ny, f_M, f_ext+idP*Nfield, ftype, f_P);
                    break;
            }

            dg_real vis_M[Nfield], vis_P[Nfield];

            vis_func(f_M, vis_M);
            vis_func(f_P, vis_P);

            const int sk = m*Nfield;
            for(fld=0;fld<Nfield;fld++){
                dg_real vis_delta = (vis_M[fld] - vis_P[fld])*0.5;

                pxFlux_M[sk+fld] = -fsc*nx*(  vis_delta );// + DELTA_2*pn_delta );
                pxFlux_P[sk+fld] =  fsc*nx*( -vis_delta );// + DELTA_2*pn_delta );

                pyFlux_M[sk+fld] = -fsc*ny*(  vis_delta );// + DELTA_2*pn_delta );
                pyFlux_P[sk+fld] =  fsc*ny*( -vis_delta );// + DELTA_2*pn_delta );
            }

        }

        for(n=0;n<Np;n++){
            const dg_real *ptLIFT = f_LIFT + n*Nfptotal;
            dg_real *px_rhsM = px + Nfield*(n+k1*Np);
            dg_real *py_rhsM = py + Nfield*(n+k1*Np);

            for(m=0;m<Nfp;m++){
                const int col1 = fp_M[m];
                const dg_real L = ptLIFT[col1];
                const dg_real *tx = pxFlux_M+m*Nfield;
                const dg_real *ty = pyFlux_M+m*Nfield;
                for(fld=0;fld<Nfield;fld++) {
                    px_rhsM[fld] += L*tx[fld];
                    py_rhsM[fld] += L*ty[fld];
                }
            }
            if (ftype == FACE_INNER){
                dg_real *px_rhsP = px + Nfield*(n+k1*Np);
                dg_real *py_rhsP = py + Nfield*(n+k1*Np);
                for(m=0;m<Nfp;m++){
                    const int col2 = fp_P[m];
                    const dg_real L = ptLIFT[col2];
                    const dg_real *sx = pxFlux_P+m*Nfield;
                    const dg_real *sy = pyFlux_P+m*Nfield;
                    for(fld=0;fld<Nfield;fld++) {
                        px_rhsP[fld] += L*sx[fld];
                        py_rhsP[fld] += L*sy[fld];
                    }
                }
            }
        }
    }
    return;
}

static void ldg_phys_vol_opt2d(dg_phys *phys, dg_phys_LDG *ldg){

    dg_cell *cell = dg_phys_cell(phys);
    dg_grid *grid = dg_phys_grid(phys);
    dg_region *region = dg_phys_region(phys);

    const int K = dg_grid_K(grid);
    const int Np = dg_cell_Np(cell);
    const int Nfield = dg_phys_Nfield(phys);

    dg_real *f_Dr = dg_cell_f_Dr(cell); //phys->cell->f_Dr;
    dg_real *f_Ds = dg_cell_f_Ds(cell); //phys->cell->f_Ds;
    dg_real *f_rhsQ = dg_phys_f_rhsQ(phys); //phys->f_rhsQ;
    dg_real *px = dg_phys_ldg_px(ldg);
    dg_real *py = dg_phys_ldg_py(ldg);

    register unsigned int k,n,m,fld;
    dg_real **drdx_p = region->drdx;
    dg_real **drdy_p = region->drdy;
    dg_real **dsdx_p = region->dsdx;
    dg_real **dsdy_p = region->dsdy;

    dg_real *rhs = vector_real_create(Nfield);
    int rhsid = 0;
    for(k=0;k<K;k++){
        // calculate viscosity term

        for(n=0;n<Np;n++){
            const dg_real *ptDr = f_Dr+n*Np; // n-th row of Dr
            const dg_real *ptDs = f_Ds+n*Np; // n-th row of Ds

            const dg_real drdx = drdx_p[k][n]; // volume geometry for n-th point
            const dg_real drdy = drdy_p[k][n]; // volume geometry for n-th point
            const dg_real dsdx = dsdx_p[k][n]; // volume geometry for n-th point
            const dg_real dsdy = dsdy_p[k][n]; // volume geometry for n-th point

            // initialize rhs
            for(m=0;m<Nfield;m++){ rhs[m] = 0; }

            for(m=0;m<Np;++m){
                const dg_real dr = ptDr[m]; // m-th column for Dr
                const dg_real ds = ptDs[m]; // m-th column for Ds
                const dg_real dx = drdx*dr+dsdx*ds;
                const dg_real dy = drdy*dr+dsdy*ds;

                const int sk = (k*Np+m)*Nfield;

                dg_real *pxk = px + sk; // variable in k-th element
                dg_real *pyk = py + sk; // variable in k-th element

                for(fld=0;fld<Nfield;fld++){ rhs[fld] += dx * pxk[fld] + dy * pyk[fld]; }
            }
            for(fld=0;fld<Nfield;fld++){ f_rhsQ[ rhsid++ ] += rhs[fld]; }
        }
    }
    vector_real_free(rhs);
    return;
}

static void ldg_phys_sruf_opt2d(dg_phys *phys, dg_phys_LDG *ldg){

    dg_edge *edge = dg_phys_edge(phys);
    dg_cell *cell = dg_phys_cell(phys);
    const int Nfield = dg_phys_Nfield(phys);
    const int Nedge = dg_edge_Nedge(edge);
    const int Nfptotal = dg_cell_Nfptotal(cell);
    const int Np = dg_cell_Np(cell);

//    dg_real *f_Q = dg_phys_f_Q(phys);
    dg_real *f_LIFT  = dg_cell_f_LIFT(cell); //phys->cell->f_LIFT;
    dg_real *f_rhsQ = dg_phys_f_rhsQ(phys);
//    dg_real *f_inQ = dg_phys_f_recvQ(phys); // phys->f_recvQ;
//    dg_real *f_ext = dg_phys_f_extQ(phys); // phys->f_extQ;

    dg_real *px = dg_phys_ldg_px(ldg);
    dg_real *py = dg_phys_ldg_py(ldg);
    dg_real *px_recv = dg_phys_ldg_x_recv(ldg);
    dg_real *py_recv = dg_phys_ldg_y_recv(ldg);

    register int f,m,n,fld,surfid=0,nodeid=0;

    for(f=0;f<Nedge;f++){
        const int k1 = edge->surfinfo[surfid++];
        const int k2 = edge->surfinfo[surfid++];
        const int f1 = edge->surfinfo[surfid++];
        const int f2 = edge->surfinfo[surfid++];
        const int ftype = edge->surfinfo[surfid++];
        const int Nfp = dg_cell_Nfp(cell)[f1];

        dg_real flux_M[Nfp*Nfield], flux_P[Nfp*Nfield];
        int fp_M[Nfp], fp_P[Nfp];
        for(m=0;m<Nfp;m++){
            const int idM = (int)edge->nodeinfo[nodeid++];
            const int idP = (int)edge->nodeinfo[nodeid++];
            fp_M[m] = (int)edge->nodeinfo[nodeid++];
            fp_P[m] = (int)edge->nodeinfo[nodeid++];
            const dg_real nx = edge->nodeinfo[nodeid++];
            const dg_real ny = edge->nodeinfo[nodeid++];
            const dg_real fsc = edge->nodeinfo[nodeid++];

            dg_real px_M[Nfield], px_P[Nfield];
            dg_real py_M[Nfield], py_P[Nfield];
            // local face2d values
            for(fld=0;fld<Nfield;fld++) {
                px_M[fld] = px[idM*Nfield+fld];
                py_M[fld] = py[idM*Nfield+fld];
            }

            // adjacent nodal values
            switch (ftype){
                case FACE_INNER:
                    for(fld=0;fld<Nfield;fld++){
                        px_P[fld] = px[idP*Nfield+fld];
                        py_P[fld] = py[idP*Nfield+fld];
                    }
                    break;
                case FACE_PARALL:
                    for(fld=0;fld<Nfield;fld++){
                        px_P[fld] = px_recv[idP*Nfield+fld];
                        py_P[fld] = py_recv[idP*Nfield+fld];
                    }
                    break;
                default: /// open boundary condition
                    for(fld=0;fld<Nfield;fld++){
                        px_P[fld] = px_M[fld];
                        py_P[fld] = py_M[fld];
                    }
                    break;
            }

            const int sk = m*Nfield;
            for(fld=0;fld<Nfield;fld++){
                dg_real px_delta = (px_M[fld] - px_P[fld])*0.5;
                dg_real py_delta = (py_M[fld] - py_P[fld])*0.5;
                dg_real pn_delta = nx*px_delta + ny*py_delta;

                flux_M[sk+fld] = -fsc*( pn_delta );
                flux_P[sk+fld] = -fsc*( pn_delta );
            }

        }

        for(n=0;n<Np;n++){
            const dg_real *ptLIFT = f_LIFT + n*Nfptotal;
            dg_real *f_rhsM = f_rhsQ + Nfield*(n+k1*Np);

            for(m=0;m<Nfp;m++){
                const int col1 = fp_M[m];
                const dg_real L = ptLIFT[col1];
                const dg_real *t = flux_M+m*Nfield;
                for(fld=0;fld<Nfield;fld++) {f_rhsM[fld] += L*t[fld];}
            }
            if (ftype == FACE_INNER){
                dg_real *f_rhsP = f_rhsQ + Nfield*(n+k2*Np);
                for(m=0;m<Nfp;m++){
                    const int col2 = fp_P[m];
                    const dg_real L = ptLIFT[col2];
                    const dg_real *s = flux_P+m*Nfield;
                    for(fld=0;fld<Nfield;fld++) {f_rhsP[fld] += L*s[fld];}
                }
            }
        }
    }
    return;
}