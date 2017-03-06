/**
 * @brief
 * Two dimensional convection problem
 *
 * @details
 * Two dimensional scalar convection problem
 * \f[\frac{\partial C}{\partial t}+\frac{\partial uC}{\partial x}+\frac{\partial vC}{\partial y}=0\f]
 *
 * Usages:
 * Use the 2 order basis with uniform triangular mesh (0) of 80 elements on each edge:
 *
 *     mpirun -n 2 -host localhost ./convection2d
 *
 * @author
 * Li Longxiang, Tianjin University, li12242@tju.edu.cn
 */
#define DEBUG 0
#include "conv_driver2d.h"
#include "conv_getparameter.h"
#include "conv_intilization.h"
#include "conv_mesh.h"
#include "conv_output.h"
#include "MultiRegions/mr_mesh_bc.h"
#include "conv_run.h"
#include "conv_extsol.h"

conv_solver2d solver; ///< global solver

/* finalize */
void conv_finalize(physField *phys);

int main(int argc, char **argv){

    /* initialize MPI */
    MPI_Init(&argc, &argv);
#if DEBUG
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    if(!procid) printf("step into conv_getparameter\n");
#endif
    /* get parameters */
    conv_getparameter(argc, argv);

    /* physical field */
    dg_cell *cell = dg_cell_creat(solver.N, solver.celltype);
#if DEBUG
    if(!procid) printf("step into conv_mesh\n");
#endif
    parallMesh *mesh = conv_mesh(cell);
    physField *phys = pf_create(3, mesh);

#if DEBUG
    if(!procid) printf("step out conv_intilization\n");
#endif
    /* initialization */
    conv_intilization(phys);

    /* set output */
    conv_setoutput(phys);

    /* run solver */
    conv_run(phys);

    /* cal norm error */
    conv_normerr(phys);

    /* finalize */
    conv_finalize(phys);

    return 0;
}

void conv_finalize(physField *phys){

    extern conv_solver2d solver;
    nc_file_close(solver.outfile);
    nc_file_free(solver.outfile);

    mr_grid_free(phys->grid);
    mr_reg_free(phys->region);
    mr_mesh_del_bc2d(phys->mesh);
    mr_mesh_free(phys->mesh);
    pf_free(phys);
    MPI_Finalize();
}