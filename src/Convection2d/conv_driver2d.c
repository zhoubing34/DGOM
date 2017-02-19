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

    /* get parameters */
    conv_getparameter(argc, argv);

    /* physical field */
    stdCell *cell = sc_create(solver.N, solver.celltype);
    parallMesh *mesh = conv_mesh(cell);
    physField *phys = pf_create(3, mesh);

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
    mr_mesh_deleteBoundary2d(phys->mesh);
    mr_mesh_free(phys->mesh);
    pf_free(phys);
    MPI_Finalize();
}