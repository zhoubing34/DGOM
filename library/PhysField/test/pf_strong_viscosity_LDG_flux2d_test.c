//
// Created by li12242 on 17/1/5.
//

#include <MultiRegions/Grid/dg_grid.h>
#include <MultiRegions/Mesh/mr_mesh.h>
#include "pf_strong_viscosity_LDG_flux2d_test.h"
#include "pf_test.h"
#include "PhysField/pf_add_LDG_solver.h"
#include "PhysField/pf_fetchBuffer.h"

static int wall_func(dg_real nx, dg_real ny, dg_real *varM, dg_real *varP,
              dg_real *pxM, dg_real *pxP, dg_real *pyM, dg_real *pyP){
    varP[0] = varM[0];
    varP[1] = varM[1];
    pxP[0] = pxM[0];
    pxP[1] = pxM[1];
    pyP[0] = pyM[0];
    pyP[1] = pyM[1];
    return 0;
}

int phys_strong_viscosity_LDG_flux2d_test(physField *phys, int verbose){
    int fail = 0;
    extern int Nfield;

    dg_mesh *mesh = phys->mesh;
    dg_region *region = phys->region;
    dg_cell *shape = phys->cell;

    const int K = phys->grid->K;
    const int Np = shape->Np;
    const int nprocs = phys->mesh->nprocs;

    int k,i;
    dg_real miu = 0.5;
    dg_real px_ext[Np*Nfield*K];
    dg_real py_ext[Np*Nfield*K];
    dg_real rhs_ext[Np*Nfield*K];

    pf_add_LDG_solver(phys);

    // assignment
    int sk = 0;
    for(k=0;k<K;k++){
        for(i=0;i<Np;i++){
            dg_real u = region->x[k][i];
            dg_real v = region->y[k][i];

            phys->f_rhsQ[sk] = 0;
            phys->f_Q[sk] = u*u; //
            phys->viscosity->vis_Q[sk] = miu;
            px_ext[sk] = 2*u*miu;
            py_ext[sk] = 0.0*miu;
            rhs_ext[sk++] = 2.0*miu*2;

            phys->f_rhsQ[sk] = 0;
            phys->f_Q[sk] = u*v+v*v; // field 1 of -(dEdx + dGdy)
            phys->viscosity->vis_Q[sk] = miu;
            px_ext[sk] = v*miu;
            py_ext[sk] = (u+2*v)*miu;
            rhs_ext[sk++] = 2.0*miu*2;
        }
    }

    MPI_Request mpi_out_requests[nprocs];
    MPI_Request mpi_in_requests[nprocs];
    int Nmess;
    pf_fetchNodeBuffer2d(phys, mpi_out_requests, mpi_in_requests, &Nmess);
    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_in_requests, instatus);

    double clockT1 = MPI_Wtime();
    pf_strong_viscosity_LDG_flux2d(phys, wall_func, wall_func, 0, 0, 0);
    double clockT2 = MPI_Wtime();
    pf_strong_viscosity_LDG_flux2d(phys, wall_func, wall_func, 0, 0, 0);

    if(!phys->mesh->procid){
        vector_double_test(__FUNCTION__, phys->viscosity->px_Q, px_ext, Np * Nfield * K);
        vector_double_test(__FUNCTION__, phys->viscosity->py_Q, py_ext, Np * Nfield * K);
        vector_double_test(__FUNCTION__, phys->f_rhsQ, rhs_ext, Np * Nfield * K);
    }

    if(verbose){
        FILE *fp = create_log(__FUNCTION__, mesh->procid, mesh->nprocs);
        fprintf(fp, "K = %d\n", phys->grid->K);
        fprintf(fp, "Nfield = %d\n", phys->Nfield);
        fprintf(fp, "Np = %d\n", phys->cell->Np);
        print_int_matrix2file(fp, "EToV", phys->grid->EToV, K, phys->cell->Nv);
        print_double_vector2file(fp, "vx", phys->grid->vx, phys->grid->Nv);
        print_double_vector2file(fp, "vy", phys->grid->vy, phys->grid->Nv);
        print_double_vector2file(fp, "f_Q", phys->f_Q, Nfield * Np * k);
        print_double_vector2file(fp, "px_Q", phys->viscosity->px_Q, Nfield * Np * K);
        print_double_vector2file(fp, "py_Q", phys->viscosity->py_Q, Nfield * Np * K);
        print_double_vector2file(fp, "f_rhsQ", phys->f_rhsQ, Nfield * Np * k);
        fclose(fp);
    }

    pf_delete_LDG_solver(phys);

    const int procid = region->procid;
    if(!procid) {
        if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    }
    return fail;
}