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
 *     mpirun -n 2 -host localhost ./SWE2D 2 80 60 tri ParabolicBowl
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 */

#include "SWEDriver2d.h"

/* private function */
void SWEFinalize(MultiReg2d *mesh, PhysDomain2d *phys, NcFile *file);

int main(int argc, char **argv){
    /* initialize MPI */
    MPI_Init(&argc, &argv);

    /* allocate solver */
    SWE_Solver2d *solver = (SWE_Solver2d *)malloc(sizeof(SWE_Solver2d));
    /* allocate mesh */
    MultiReg2d   *mesh = SWE_Mesh2d(argv, solver);
    /* allocate physical domain */
    PhysDomain2d *phys = SWE_Init2d(argv, solver, mesh);

    /* set output file */
    NcFile *outfile = SWE_SetNcOutput2d(phys, solver);
    /* solve  */
    SWERun2d(phys, solver, outfile);
    /* finalize */
    CloseNcFile(outfile);
    SWEFinalize(mesh, phys, outfile);
    MPI_Finalize();

    return 0;
}

void SWEFinalize(MultiReg2d *mesh,
                 PhysDomain2d *phys, NcFile *file){
    FreeStdRegions2d(mesh->stdcell);
    FreeMultiReg2d(mesh);
    FreePhysDomain2d(phys);
    FreeNcFile(file);
}