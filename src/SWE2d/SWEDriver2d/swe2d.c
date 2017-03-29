/**
 * @file
 * SWEDriver2d.c
 *
 * @brief
 * Two dimensional shallow water equations.
 *
 * @details
 * Two dimensional shallow water equations
 * \f[ \begin{equation} \begin{array}{c}
 * \frac{\partial U}{\partial t} + \frac{\partial E(U)}{\partial x} + \frac{\partial G(U)}{\partial y}  = S(U) \cr
 * U = \begin{bmatrix}\eta \cr q_x \cr q_y \end{bmatrix} \quad E = \begin{bmatrix}q_x \cr
 * \frac{q_x^2}{h} + \frac{1}{2}gh^2 \cr
 * \frac{q_xq_y}{h}\end{bmatrix} \quad G = \begin{bmatrix}q_y \cr \frac{q_xq_y}{h} \cr
 * \frac{q_y^2}{h} + \frac{1}{2}gh^2 \end{bmatrix} \quad S = \begin{bmatrix}0 \cr -gh\frac{\partial z}{\partial x} \cr
 * -gh\frac{\partial z}{\partial y} \end{bmatrix}
 * \end{array} \end{equation} \f]
 *
 * Usages:
 * Use the 2 order basis with an uniform mesh of 80 and 60 elements along x and y coordinate respectively.
 * For the ParabolicBowl test case:
 *
 *     mpirun -n 2 -host localhost ./SWE2D
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 */

#include "swe2d.h"
#include "../SWELib/swe_lib.h"

static void swe_finalize();

SWE_Solver solver;

int main(int argc, char **argv){

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    /* read input file */
    swe_input(argc, argv);

    /* initialize */
    swe_init();

    /* output */
    swe_output();

    /* run */
    swe_run();

    /* finalize */
    swe_finalize();

    return 0;
}

/* finalize the SWE and deallocate the variable */
static void swe_finalize(){

    extern SWE_Solver solver;
    /* physcal field */
    dg_cell_free(dg_phys_cell(solver.phys));
    dg_grid_free(dg_phys_grid(solver.phys));
    dg_region_free(dg_phys_region(solver.phys));
    dg_mesh_free(dg_phys_mesh(solver.phys));
    dg_edge_free(dg_phys_edge(solver.phys));
    dg_phys_free(solver.phys);

    /* output file */
    nc_file_close(solver.outfile);
    nc_file_free(solver.outfile);

    /* solver */
    vector_real_free(solver.m);
    MPI_Finalize();
}