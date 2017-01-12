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

#include "swe_dirver2d.h"
#include "MultiRegions/mr_mesh_addBoundary.h"

swe_solver solver;

int main(int argc, char **argv){

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    /* get parameter */


    /* allocate region domain */
    MultiReg2d   *region = SWE_Mesh2d(argv, solver);
    /* allocate physical domain */
    PhysDomain2d *phys = SWE_Init2d(argv, solver, region);
    /* set output file */
    NcFile *outfile = SWE_SetNcOutput2d(phys, solver);
    /* solve  */
    SWE_Run2d(phys, solver, outfile);

    /* finalize */
    SWE_Finalize2d(region, phys, outfile);
    MPI_Finalize();

    return 0;
}

/* finalize the SWE and deallocate the variable */
static void swe_finalize(physField *phys){
    nc_file_close(solver.outfile);
    nc_file_free(solver.outfile);

    mr_grid_free(phys->grid);
    mr_reg_free(phys->region);
    mr_mesh_deleteBoundary2d(phys->mesh);
    mr_mesh_free(phys->mesh);
    phys_free(phys);

    MPI_Finalize();
}