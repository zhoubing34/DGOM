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

#include "conv_driver.h"

Conv_Solver solver; ///< global solver

/* finalize */
void conv_finalize();

int main(int argc, char **argv){

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    /* get parameters */
    Conv_Case_Type case_type = conv_input(argc, argv);

    /* initialization */
    conv_init(case_type);

    /* set output */
    conv_setoutput();

    /* run solver */
    conv_run();

    /* cal norm error */
    conv_normerr(case_type);

    /* finalize */
    conv_finalize();

    return 0;
}

void conv_finalize(){
    extern Conv_Solver solver;
    nc_file_close(solver.outfile);
    nc_file_free(solver.outfile);

    dg_phys *phys = solver.phys;
    dg_cell_free(phys->cell);
    dg_grid_free(phys->grid);
    dg_region_free(phys->region);
    dg_mesh_free(phys->mesh);
    dg_phys_free(phys);
    MPI_Finalize();
}