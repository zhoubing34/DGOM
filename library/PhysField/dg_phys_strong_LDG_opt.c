//
// Created by li12242 on 17/1/1.
//

#include "dg_phys_strong_LDG_opt.h"

#define DEBUG 0

#if DEBUG
#include "LibUtilities/UTest.h"
#endif

static void ldg_auxi_vol_opt2d(dg_phys *phys, dg_phys_LDG *ldg, Vis_Fun vis_func);
static void ldg_auxi_surf_opt2d(dg_phys *phys, dg_phys_LDG *ldg, Vis_Fun vis_func,
                                Wall_Condition slipwall_func, Wall_Condition nonslipwall_func, OBC_Fun obc_fun);
static void ldg_phys_vol_opt2d(dg_phys *phys, dg_phys_LDG *ldg);

void dg_phys_LDG_solve_vis_opt2d(dg_phys *phys, Vis_Fun vis_func,
                                 Wall_Condition slipwall_func,
                                 Wall_Condition non_slipwall_func,
                                 OBC_Fun obc_fun){

    dg_mesh *mesh = dg_phys_mesh(phys);
    dg_phys_LDG *ldg = phys->ldg;
    const int nprocs = dg_mesh_nprocs(mesh);

    /* 0. send/recv the viscosity value */
    const int Nfield = dg_phys_Nfield(phys);
    int Nmess;
    MPI_Request mpi_x_send_requests[nprocs], mpi_x_recv_requests[nprocs];
    MPI_Request mpi_y_send_requests[nprocs], mpi_y_recv_requests[nprocs];
    Nmess = mesh->fetch_node_buffer(mesh, Nfield, ldg->sqrt_miux, ldg->px_recv,
                                    mpi_x_send_requests, mpi_x_recv_requests);
    Nmess = mesh->fetch_node_buffer(mesh, Nfield, ldg->sqrt_miuy, ldg->py_recv,
                                    mpi_y_send_requests, mpi_y_recv_requests);

    /* 1. calculate the auxiliary variables - volume integral */
    ldg_auxi_vol_opt2d(phys, ldg, vis_func);
    /* wait for message passing */
    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_x_send_requests, instatus);
    MPI_Waitall(Nmess, mpi_x_recv_requests, instatus);
    MPI_Waitall(Nmess, mpi_y_send_requests, instatus);
    MPI_Waitall(Nmess, mpi_y_recv_requests, instatus);
    /* 2. calculate the auxiliary variables - surface integral */
    ldg_auxi_surf_opt2d(phys, ldg, vis_func, slipwall_func, non_slipwall_func, obc_fun);

    /* 3. calculate the volume term */
    ldg_phys_vol_opt2d(phys, ldg);
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
    dg_real *miu_xp = dg_phys_ldg_miux(ldg);
    dg_real *miu_yp = dg_phys_ldg_miuy(ldg);

    dg_real *vis_term = vector_real_create(Np*Nfield);
    dg_real *px_rhs = vector_real_create(Nfield);
    dg_real *py_rhs = vector_real_create(Nfield);

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

            const dg_real *miu_x = miu_xp + (k*Np + n)*Nfield;
            const dg_real *miu_y = miu_yp + (k*Np + n)*Nfield;

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
                    px_rhs[fld] += miu_x[fld]*dx*vis[fld];
                    py_rhs[fld] += miu_y[fld]*dy*vis[fld];
                }
            }

            const int sk = (k*Np+n)*Nfield;
            for(fld=0;fld<Nfield;fld++){
                px[ sk+fld ] = px_rhs[fld];
                py[ sk+fld ] = py_rhs[fld];
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
    dg_real *miu_xp = dg_phys_ldg_miux(ldg);
    dg_real *miu_yp = dg_phys_ldg_miuy(ldg);
    dg_real *miux_recv = dg_phys_ldg_x_recv(ldg);
    dg_real *miuy_recv = dg_phys_ldg_y_recv(ldg);

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
            dg_real miux_M[Nfield], miux_P[Nfield];
            dg_real miuy_M[Nfield], miuy_P[Nfield];
            // local face2d values
            for(fld=0;fld<Nfield;fld++) {
                f_M[fld] = f_Q[idM*Nfield+fld];
                miux_M[fld] = miu_xp[idM*Nfield+fld];
                miuy_M[fld] = miu_yp[idM*Nfield+fld];
            }

            // default adjacent node viscosity value
            for(fld=0;fld<Nfield;fld++){
                miux_P[fld] = miu_xp[idM*Nfield+fld];
                miuy_P[fld] = miu_yp[idM*Nfield+fld];
            }
            // adjacent nodal values
            switch (ftype){
                case FACE_INNER:
                    for(fld=0;fld<Nfield;fld++){
                        f_P[fld] = f_Q[idP*Nfield+fld];
                        miux_P[fld] = miu_xp[idP*Nfield+fld];
                        miuy_P[fld] = miu_yp[idP*Nfield+fld];
                    }
                    break;
                case FACE_PARALL:
                    for(fld=0;fld<Nfield;fld++){
                        f_P[fld] = f_inQ[idP*Nfield+fld];
                        miux_P[fld] = miux_recv[idP*Nfield+fld];
                        miuy_P[fld] = miuy_recv[idP*Nfield+fld];
                    }
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
                dg_real miux_mean = (miux_M[fld]+miux_P[fld])*0.5;
                dg_real miuy_mean = (miuy_M[fld]+miuy_P[fld])*0.5;

                pxFlux_M[sk+fld] = -fsc*nx*(miux_M[fld]*vis_M[fld] - miux_mean*vis_P[fld] );
                pxFlux_P[sk+fld] = +fsc*nx*(miux_P[fld]*vis_P[fld] - miux_mean*vis_M[fld] );

                pyFlux_M[sk+fld] = -fsc*ny*(miuy_M[fld]*vis_M[fld] - miuy_mean*vis_P[fld] );
                pyFlux_P[sk+fld] = +fsc*ny*(miuy_P[fld]*vis_P[fld] - miuy_mean*vis_M[fld] );
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
    dg_real *miu_xp = dg_phys_ldg_miux(ldg);
    dg_real *miu_yp = dg_phys_ldg_miuy(ldg);

    dg_real *rhs = vector_real_create(Nfield);
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

                const dg_real *miu_x = miu_xp + (k*Np + m)*Nfield;
                const dg_real *miu_y = miu_yp + (k*Np + m)*Nfield;

                const int sk = (k*Np+m)*Nfield;
                dg_real *pxk = px + sk; // variable in k-th element
                dg_real *pyk = py + sk; // variable in k-th element

                for(fld=0;fld<Nfield;fld++){
                    rhs[fld] += dx * miu_x[fld] * pxk[fld]
                            + dy * miu_y[fld] * pyk[fld];
                }
            }
            const int sk = (k*Np+n)*Nfield;
            for(fld=0;fld<Nfield;fld++){ f_rhsQ[ sk+fld ] += rhs[fld]; }
        }
    }
    vector_real_free(rhs);
    return;
}

//static void ldg_phys_surf_opt2d(dg_phys *phys, dg_phys_LDG *ldg){
//
//    return;
//}


//void pf_strong_viscosity_LDG_flux2d(physField *phys,
//                                    viscosity_wall_condition_func slipwall_condition,
//                                    viscosity_wall_condition_func non_slipwall_condition,
//                                    dg_real c11, dg_real c12, dg_real c22){
//
//    // check LDG solver is initialized
//    if(phys->viscosity == NULL){
//        fprintf(stderr, "PhysField (%s): line %d\n"
//                "The LDG solver for viscosity is not allocated\n", __FILE__, __LINE__);
//    }
//
//    const int procid = phys->grid->procid;
//    const int nprocs = phys->grid->nprocs;
//    const int Nfp = phys->cell->Nfp;
//    const int Nfield = phys->Nfield;
//
//    dg_real *p_Q = phys->viscosity->px_Q;
//    dg_real *q_Q = phys->viscosity->py_Q;
//    dg_real *p_inQ = phys->viscosity->px_inQ;
//    dg_real *q_inQ = phys->viscosity->py_inQ;
//    dg_real *p_outQ = phys->viscosity->px_outQ;
//    dg_real *q_outQ = phys->viscosity->py_outQ;
//
//    register int n;
//
//    // calculate the auxiliary variable
//    phys_auxiliaryflux2d(phys, slipwall_condition, non_slipwall_condition, c11, c12, c22);
//
//
//
//    // fetch p_Q and q_Q buffer
//    for(n=0;n<phys->parallNodeNum;++n){
//        p_outQ[n] = p_Q[phys->nodeIndexOut[n]];
//        q_outQ[n] = q_Q[phys->nodeIndexOut[n]];
//    }
//
//    int Nout[nprocs];
//    for(n=0;n<nprocs;n++){
//        Nout[n] = phys->mesh->Parf[n]*Nfield*Nfp;
//    }
//
//    MPI_Request mpi_send_requests[nprocs], mpi_recv_requests[nprocs];
//    int Nmess=0;
//
//    /* do sends and recv for p_Q */
//    pf_fetchBuffer(procid, nprocs, Nout, p_outQ, p_inQ,
//                   mpi_send_requests, mpi_recv_requests, &Nmess);
//
//    MPI_Status instatus[nprocs];
//    MPI_Waitall(Nmess, mpi_recv_requests, instatus);
//    MPI_Waitall(Nmess, mpi_send_requests, instatus);
//
//    /* do sends and recv for q_Q */
//    pf_fetchBuffer(procid, nprocs, Nout, q_outQ, q_inQ,
//                   mpi_send_requests, mpi_recv_requests, &Nmess);
//    MPI_Waitall(Nmess, mpi_recv_requests, instatus);
//    MPI_Waitall(Nmess, mpi_send_requests, instatus);
//
//    // calculate the RHS of physical field
//    phys_viscosityflux2d(phys, slipwall_condition, non_slipwall_condition, c11, c12, c22);
//
//}
//
//
//static void phys_auxiliaryflux2d(physField *phys,
//                                 viscosity_wall_condition_func slipwall_condition,
//                                 viscosity_wall_condition_func non_slipwall_condition,
//                                 dg_real c11, dg_real c12, dg_real c22){
//
//    const int K = phys->grid->K;
//    const int Np = phys->cell->Np;
//    const int Nfp = phys->cell->Nfp;
//    const int Nfaces = phys->cell->Nfaces;
//    const int Nfield = phys->Nfield;
//
//    dg_real *p_Q = phys->viscosity->px_Q;
//    dg_real *q_Q = phys->viscosity->py_Q;
//    dg_real *p_inQ = phys->viscosity->px_inQ;
//    dg_real *q_inQ = phys->viscosity->py_inQ;
//
//    dg_real *f_inQ = phys->f_inQ;
//    dg_real *f_ext = phys->f_ext;
//
//    dg_real *f_Q = phys->f_Q;
//    dg_real *f_Dr = phys->cell->f_Dr;
//    dg_real *f_Ds = phys->cell->f_Ds;
//    dg_real *f_LIFT = phys->cell->f_LIFT;
//    dg_real *vgeo = phys->vgeo;
//    dg_real *surfinfo = phys->surfinfo;
//
//    dg_real *vis_Q = phys->viscosity->vis_Q;
//
//    register int m,n,k,fld,geoid=0,surfid=0;
//
//    dg_real f_varM[Nfield], f_varP[Nfield];
//    dg_real p_varM[Nfield], p_varP[Nfield], p_dflux[Nfp*Nfaces*Nfield];
//    dg_real q_varM[Nfield], q_varP[Nfield], q_dflux[Nfp*Nfaces*Nfield];
//    dg_real vis_varM[Nfield];
//
//#if DEBUG
//    FILE *fp = CreateLog("phys_auxiliaryflux2d_log",
//                         phys->mesh->procid, phys->mesh->nprocs);
//#endif
//
//    for(k=0;k<K;k++){
//        // volume integral for p and q
//        for(n=0;n<Np;n++){
//            const dg_real *ptDr = f_Dr+n*Np; // n-th row of Dr
//            const dg_real *ptDs = f_Ds+n*Np; // n-th row of Ds
//
//            const dg_real drdx = vgeo[geoid++]; // volume geometry for n-th point
//            const dg_real drdy = vgeo[geoid++]; // volume geometry for n-th point
//            const dg_real dsdx = vgeo[geoid++]; // volume geometry for n-th point
//            const dg_real dsdy = vgeo[geoid++]; // volume geometry for n-th point
//            const dg_real *visc = phys->viscosity->vis_Q+(k*Np+n)*Nfield;
//
//            dg_real *p = p_Q + (n+k*Np)*Nfield;
//            dg_real *q = q_Q + (n+k*Np)*Nfield;
//
//            // initialization: set px and py to zero
//            for(fld=0;fld<Nfield;fld++){
//                p[fld] = 0;
//                q[fld] = 0;
//            }
//
//            for(m=0;m<Np;++m){
//                const dg_real dr = ptDr[m]; // m-th column for Dr
//                const dg_real ds = ptDs[m]; // m-th column for Ds
//                const dg_real dx = drdx*dr+dsdx*ds;
//                const dg_real dy = drdy*dr+dsdy*ds;
//
//                const dg_real *c = f_Q + (m+k*Np)*Nfield;
//
//                for(fld=0;fld<Nfield;fld++){
//                    p[fld] += visc[fld]*dx*c[fld];
//                    q[fld] += visc[fld]*dy*c[fld];
//                }
//            }
//        }
//
//        // surface integral
//        for(n=0;n<Nfp*Nfaces;n++){
//            int  idM = (int)surfinfo[surfid++];
//            int  idP = (int)surfinfo[surfid++];
//            const dg_real fsc = surfinfo[surfid++];
//            const int  bstype = (int)surfinfo[surfid++];
//            const dg_real nx = surfinfo[surfid++];
//            const dg_real ny = surfinfo[surfid++];
//
//            // local face2d values
//            for(fld=0;fld<Nfield;fld++){
//                f_varM[fld] = f_Q[idM];
//                p_varM[fld] = p_Q[idM];
//                q_varM[fld] = q_Q[idM];
//                vis_varM[fld] = vis_Q[idM++];
//            }
//
//#if DEBUG
//            fprintf(fp, "\nNfp=%d, bctype=%d, ", n, bstype);
//#endif
//
//            // adjacent nodal values
//            switch (bstype){
//                case INNERLOC:
//                    for(fld=0;fld<Nfield;fld++){
//                        f_varP[fld] = f_Q[idP];
//                        p_varP[fld] = p_Q[idP];
//                        q_varP[fld] = q_Q[idP++];
//                    }
//                    break;
//                case INNERBS:
//                    for(fld=0;fld<Nfield;fld++){
//                        f_varP[fld] = f_inQ[idP];
//                        p_varP[fld] = p_inQ[idP];
//                        q_varP[fld] = q_inQ[idP++];
//                    }
//                    break;
//                case SLIPWALL:
//                    slipwall_condition(nx, ny, f_varM, f_varP,
//                                       p_varM, p_varP, q_varM, q_varP);
//                    break;
//                case NSLIPWALL:
//                    non_slipwall_condition(nx, ny, f_varM, f_varP,
//                                           p_varM, p_varP, q_varM, q_varP);
//                    break;
//                default: // open boundary condition
//                    for(fld=0;fld<Nfield;fld++){
//                        f_varP[fld] = f_ext[idP++];
//                        p_varP[fld] = p_varM[fld]; // zero gradient
//                        q_varP[fld] = q_varM[fld];
//                    }
//                    break;
//            }
//
//            dg_real *flux_P = p_dflux + n*Nfield;
//            dg_real *flux_Q = q_dflux + n*Nfield;
//
//            for(fld=0;fld<Nfield;fld++){
//                const dg_real df = (f_varM[fld] - f_varP[fld]);
//                const dg_real dp = nx*(p_varM[fld] - p_varP[fld]) + ny*(q_varM[fld] - q_varP[fld]);
//#if DEBUG
//                fprintf(fp, "%f, ", df);
//#endif
//                flux_P[fld] = vis_varM[fld]*fsc*( nx*( -df*0.5 + c12*nx*df - c22*dp ) );
//                flux_Q[fld] = vis_varM[fld]*fsc*( ny*( -df*0.5 + c12*ny*df - c22*dp ) );
//            }
//        }
//
//#if DEBUG
//        PrintVector2File(fp, "elemental flux_P", p_dflux, Nfield*Nfp*Nfaces);
//        PrintVector2File(fp, "elemental flux_Q", q_dflux, Nfield*Nfp*Nfaces);
//#endif
//
//        for(n=0;n<Np;n++){
//            const dg_real *ptLIFT = f_LIFT + n*Nfp*Nfaces;
//
//            dg_real *p_rhsQ = p_Q + Nfield*(n+k*Np);
//            dg_real *q_rhsQ = q_Q + Nfield*(n+k*Np);
//            for(m=0;m<Nfp*Nfaces;m++){
//                const dg_real L = ptLIFT[m];
//                dg_real *flux_P = p_dflux+m*Nfield;
//                dg_real *flux_Q = q_dflux+m*Nfield;
//
//                for(fld=0;fld<Nfield;fld++){
//                    p_rhsQ[fld] += L*flux_P[fld];
//                    q_rhsQ[fld] += L*flux_Q[fld];
//                }
//
//            }
//        }
//    }
//#if DEBUG
//    fclose(fp);
//#endif
//
//}
//
//static void phys_viscosityflux2d(physField *phys,
//                                 viscosity_wall_condition_func slipwall_condition,
//                                 viscosity_wall_condition_func non_slipwall_condition,
//                                 dg_real c11, dg_real c12, dg_real c22){
//    const int K = phys->grid->K;
//    const int Np = phys->cell->Np;
//    const int Nfp = phys->cell->Nfp;
//    const int Nfaces = phys->cell->Nfaces;
//    const int Nfield = phys->Nfield;
//
//    dg_real *f_Q = phys->f_Q;
//    dg_real *f_Dr = phys->cell->f_Dr;
//    dg_real *f_Ds = phys->cell->f_Ds;
//    dg_real *f_LIFT = phys->cell->f_LIFT;
//    dg_real *f_rhs = phys->f_rhsQ;
//    dg_real *vgeo = phys->vgeo;
//    dg_real *surfinfo = phys->surfinfo;
//
//    dg_real *p_Q = phys->viscosity->px_Q;
//    dg_real *q_Q = phys->viscosity->py_Q;
//    dg_real *p_inQ = phys->viscosity->px_inQ;
//    dg_real *q_inQ = phys->viscosity->py_inQ;
//
//    dg_real *f_inQ = phys->f_inQ;
//    dg_real *f_ext = phys->f_ext;
//
//    register int k,n,m,fld,geoid=0,surfid=0;
//
//    dg_real f_varM[Nfield], f_varP[Nfield], f_dflux[Nfp*Nfaces*Nfield];
//    dg_real p_varM[Nfield], p_varP[Nfield];
//    dg_real q_varM[Nfield], q_varP[Nfield];
//
//#if DEBUG
//    FILE *fp = CreateLog("phys_viscosityflux2d_log",
//                         phys->mesh->procid, phys->mesh->nprocs);
//#endif
//
//    for(k=0;k<K;k++){
//        // volume integral for p and q
//        for(n=0;n<Np;n++){
//            const dg_real *ptDr = f_Dr+n*Np; // n-th row of Dr
//            const dg_real *ptDs = f_Ds+n*Np; // n-th row of Ds
//
//            const dg_real drdx = vgeo[geoid++]; // volume geometry for n-th point
//            const dg_real drdy = vgeo[geoid++]; // volume geometry for n-th point
//            const dg_real dsdx = vgeo[geoid++]; // volume geometry for n-th point
//            const dg_real dsdy = vgeo[geoid++]; // volume geometry for n-th point
//
//            dg_real *rhs = f_rhs + (k*Np+n)*Nfield;
//
//            for(m=0;m<Np;++m){
//                const dg_real dr = ptDr[m]; // m-th column for Dr
//                const dg_real ds = ptDs[m]; // m-th column for Ds
//                const dg_real dx = drdx*dr+dsdx*ds;
//                const dg_real dy = drdy*dr+dsdy*ds;
//
//                const dg_real *p = p_Q + (m+k*Np)*Nfield;
//                const dg_real *q = q_Q + (m+k*Np)*Nfield;
//
//                for(fld=0;fld<Nfield;fld++){
//                    rhs[fld] += dx*p[fld] + dy*q[fld];
//                }
//            }
//        }
//
//        for(n=0;n<Nfaces*Nfp;n++){
//            int  idM = (int)surfinfo[surfid++];
//            int  idP = (int)surfinfo[surfid++];
//            const dg_real fsc = surfinfo[surfid++];
//            const int  bstype = (int)surfinfo[surfid++];
//            const dg_real nx = surfinfo[surfid++];
//            const dg_real ny = surfinfo[surfid++];
//
//            // local face2d values
//            for(fld=0;fld<Nfield;fld++){
//                f_varM[fld] = f_Q[idM];
//                p_varM[fld] = p_Q[idM];
//                q_varM[fld] = q_Q[idM++];
//            }
//
//#if DEBUG
//            fprintf(fp, "\nNfp=%d, bctype=%d, ", n, bstype);
//#endif
//
//            // adjacent nodal values
//            switch (bstype){
//                case INNERLOC:
//                    for(fld=0;fld<Nfield;fld++){
//                        f_varP[fld] = f_Q[idP];
//                        p_varP[fld] = p_Q[idP];
//                        q_varP[fld] = q_Q[idP++];
//                    }
//                    break;
//                case INNERBS:
//                    for(fld=0;fld<Nfield;fld++){
//                        f_varP[fld] = f_inQ[idP];
//                        p_varP[fld] = p_inQ[idP];
//                        q_varP[fld] = q_inQ[idP++];
//                    }
//                    break;
//                case SLIPWALL:
//                    slipwall_condition(nx, ny, f_varM, f_varP,
//                                       p_varM, p_varP, q_varM, q_varP);
//                    break;
//                case NSLIPWALL:
//                    non_slipwall_condition(nx, ny, f_varM, f_varP,
//                                           p_varM, p_varP, q_varM, q_varP);
//                    break;
//                default: // open boundary condition
//                    for(fld=0;fld<Nfield;fld++){
//                        f_varP[fld] = f_ext[idP++];
//                        p_varP[fld] = p_varM[fld]; // zero gradient
//                        q_varP[fld] = q_varM[fld];
//                    }
//                    break;
//            }
//
//            dg_real *flux = f_dflux + n*Nfield;
//
//            for(fld=0;fld<Nfield;fld++){
//                const dg_real df = (f_varM[fld] - f_varP[fld]);
//                const dg_real dp = (p_varM[fld] - p_varP[fld]);
//                const dg_real dq = (q_varM[fld] - q_varP[fld]);
//                const dg_real dpn = nx*dp + ny*dq;
//
//#if DEBUG
//                fprintf(fp, "%f, %f, %f, ", df, dp, dq);
//#endif
//                flux[fld] = fsc*( -nx*dp*0.5 - ny*dq*0.5 - c11*df - c12*(nx+ny)*dpn);
//            }
//        }
//
//#if DEBUG
//        PrintVector2File(fp, "elemental flux", f_dflux, Nfield*Nfp*Nfaces);
//#endif
//
//        for(n=0;n<Np;n++){
//            const dg_real *ptLIFT = f_LIFT + n*Nfp*Nfaces;
//            for(m=0;m<Nfp*Nfaces;m++){
//                const dg_real L = ptLIFT[m];
//                dg_real *flux_Q = f_dflux+m*Nfield;
//                dg_real *f_rhsQ = phys->f_rhsQ + Nfield*(n+k*Np);
//
//                for(fld=0;fld<Nfield;fld++)
//                    f_rhsQ[fld] += L*flux_Q[fld];
//            }
//        }
//    }
//
//#if DEBUG
//    fclose(fp);
//#endif
//}