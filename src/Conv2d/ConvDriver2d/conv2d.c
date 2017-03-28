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

#include "conv2d.h"


Conv_Solver solver; ///< global solver

/* finalize */
void conv_finalize();

int main(int argc, char **argv){

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    /* run type */
    conv_input(argc, argv);

    /* read parameter and initialize */
    conv_init();

    /* set output */
    conv_setoutput();

    /* run solver */
    conv_run();

    /* finalize */
    conv_finalize();

    return 0;
}

void conv_finalize(){
    dg_phys *phys = solver.phys;

    nc_file_close(solver.outfile);
    nc_file_free(solver.outfile);

    dg_cell_free(dg_phys_cell(phys));
    dg_grid_free(dg_phys_grid(solver.phys));
    dg_region_free(dg_phys_region(phys));
    dg_mesh_free(dg_phys_mesh(phys));
    dg_phys_free(phys);
    MPI_Finalize();
}