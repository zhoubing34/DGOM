//
// Created by li12242 on 17/1/14.
//

#include "swe_lib.h"
#define DEBUG 0
#if DEBUG
#include "Utility/unit_test.h"
#endif

static void swe_source(dg_phys *phys);

int swe_slip_wall(dg_real nx, dg_real ny, dg_real *varM, dg_real *varP){
    const dg_real hM  = varM[0];
    const dg_real qxM = varM[1];
    const dg_real qyM = varM[2];
    dg_real qnM =  qxM*nx + qyM*ny; // outward normal flux
    dg_real qvM = -qxM*ny + qyM*nx; // outward tangential flux
    // adjacent value
    varP[0] = hM;
    varP[1] = (-qnM)*nx - qvM*ny;
    varP[2] = (-qnM)*ny + qvM*nx;
    varP[3] = 0;
    return 0;
}

int swe_nonslip_wall(dg_real nx, dg_real ny, dg_real *varM, dg_real *varP){
    const dg_real hM  = varM[0];
    const dg_real qxM = varM[1];
    const dg_real qyM = varM[2];
    // adjacent value
    varP[0] = hM;
    varP[1] = -qxM;
    varP[2] = -qyM;
    varP[3] = 0;
    return 0;
}

/**
 * @brief
 * @param phys
 * @param frka
 * @param frkb
 * @param fdt
 */
void swe_rhs(dg_phys *phys, dg_real frka, dg_real frkb, dg_real fdt){

    const int K = dg_grid_K( dg_phys_grid(phys) );
    const int Np = dg_cell_Np( dg_phys_cell(phys) );
    const int Nfield = dg_phys_Nfield(phys);
    const int nprocs = dg_grid_nprocs( dg_phys_grid(phys) );

    /* mpi request buffer */
    MPI_Request *mpi_send_requests = (MPI_Request *) calloc((size_t) nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_recv_requests = (MPI_Request *) calloc((size_t) nprocs, sizeof(MPI_Request));

    int Nmess=0;
    /* fetch nodal value through all procss */
    Nmess = phys->fetch_node_buffer(phys, mpi_send_requests, mpi_recv_requests);

    /* volume integral */
    dg_phys_strong_vol_opt2d(phys, swe_flux_term);

    /* waite to recv */
    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_recv_requests, instatus);

    /* surface integral */
    dg_phys_strong_surf_opt2d(phys, swe_slip_wall, swe_nonslip_wall, swe_flux_term, swe_hll_flux);

    /* waite for finishing send buffer */
    MPI_Waitall(Nmess, mpi_send_requests, instatus);

    /* viscosity flux */
    /* source term */
    swe_source(phys);


#if DEBUG
    FILE *fp = create_log(__FUNCTION__, nprocs, dg_grid_procid(dg_phys_grid(phys)));
    print_double_vector2file(fp, "f_rhsQ", f_rhsQ, K*Np*Nfield);
    fclose(fp);
#endif
    int t,fld;
    for(t=0;t<K*Np;++t){

        dg_real *f_resQ = dg_phys_f_resQ(phys) + t*Nfield;
        dg_real *f_Q = dg_phys_f_Q(phys) + t*Nfield;
        const dg_real *f_rhsQ = dg_phys_f_rhsQ(phys) + t*Nfield;

        for(fld=0;fld<(Nfield-1);fld++){
            // calculate the resdiual of the equation
            f_resQ[fld] = frka*f_resQ[fld] + fdt*f_rhsQ[fld];
            // evaluate scalar at next internal time step
            f_Q[fld] += frkb*f_resQ[fld];
        }
    }

    free(mpi_recv_requests);
    free(mpi_send_requests);

    return;
}

/**
 * @brief add source term to the rhs field
 */
static void swe_source(dg_phys *phys){

    extern SWE_Solver solver;
    const dg_real gra = solver.gra;
    const dg_real hcrit = solver.hcrit;

    dg_cell *cell = dg_phys_cell(phys);
    dg_region *region = dg_phys_region(phys);
    const dg_real *f_Dr   = dg_cell_f_Dr(cell);
    const dg_real *f_Ds   = dg_cell_f_Ds(cell);

    const int K = dg_grid_K( dg_phys_grid(phys) );
    const int Np = dg_cell_Np(cell);
    const int Nfields  = dg_phys_Nfield(phys);

    dg_real **drdx_p = dg_region_drdx(region);
    dg_real **drdy_p = dg_region_drdy(region);
    dg_real **dsdx_p = dg_region_dsdx(region);
    dg_real **dsdy_p = dg_region_dsdy(region);

    const double p = 10/3;
    register int k,n,m;

    for(k=0;k<K;k++){
        //dg_real *bot = solver->bot + k*Np;
        dg_real *f_Q = dg_phys_f_Q(phys) + k*Np*Nfields;
        dg_real *f_rhsQ = dg_phys_f_rhsQ(phys) + k*Np*Nfields;
        for (n=0;n<Np;++n) {
            const dg_real *ptDr = f_Dr + n*Np;
            const dg_real *ptDs = f_Ds + n*Np;
            const dg_real drdx  = drdx_p[k][n], drdy = drdy_p[k][n];
            const dg_real dsdx  = dsdx_p[k][n], dsdy = dsdy_p[k][n];

            //const int ind = (k*Np+n)*Nfields;
            const dg_real h = f_Q[n*Nfields];
            const dg_real qx = f_Q[n*Nfields+1];
            const dg_real qy = f_Q[n*Nfields+2];
            const dg_real qn = sqrt(qx*qx+qy*qy);
            const dg_real Mann = solver.m[k*Np+n];
            const dg_real Mann2 = Mann*Mann;

            dg_real rhs_qx = 0, rhs_qy = 0;
            if(h>hcrit){ // for wet nodes
                for (m=0;m<Np;++m) {
                    const dg_real dr = ptDr[m];
                    const dg_real ds = ptDs[m];
                    const dg_real dx = drdx*dr + dsdx*ds;
                    const dg_real dy = drdy*dr + dsdy*ds;

                    const dg_real b = f_Q[m*Nfields+3];

                    rhs_qx += dx * b;
                    rhs_qy += dy * b;
                }

                rhs_qx += qx*qn*Mann2/pow(h, p);
                rhs_qy += qy*qn*Mann2/pow(h, p);
            }
            //const int ind = (k*Np+n)*Nfields;
            f_rhsQ[n*Nfields+1] += -gra*h*rhs_qx;
            f_rhsQ[n*Nfields+2] += -gra*h*rhs_qy;
        }
    }
    return;
}