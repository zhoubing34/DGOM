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

#include <PhysField/dg_phys.h>
#include "swe_driver2d.h"
#include "swe_init.h"
#include "swe_output.h"
#include "swe_run.h"
#include "swe_extsol.h"
static void swe_finalize(SWE_Solver*);

int main(int argc, char **argv){

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    /* initialize */
    SWE_Solver *solver = swe_init(argc, argv);

    /* output */
    solver->ncfile = swe_output(solver);

    /* run */
    swe_run(solver);

    /* norm error */
    //swe_normerr(solver);

    /* finalize */
    swe_finalize(solver);

    return 0;
}

/* finalize the SWE and deallocate the variable */
static void swe_finalize(SWE_Solver *solver){

    /* physcal field */
    dg_cell_free(solver->phys->cell);
    dg_grid_free(solver->phys->grid);
    mr_reg_free(solver->phys->region);
    mr_mesh_del_bc2d(solver->phys->mesh);
    mr_mesh_free(solver->phys->mesh);
    pf_free(solver->phys);

    /* output file */
    nc_file_close(solver->ncfile);
    nc_file_free(solver->ncfile);

    /* solver */
    vector_real_free(solver->bot);
    MPI_Finalize();
}