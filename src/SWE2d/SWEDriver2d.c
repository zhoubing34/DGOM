/**
 * @brief
 * Two dimensional shallow water equations
 *
 * @details
 * Two dimensional shallow water equations
 * \f[ \begin{equation} \begin{array}{c}
 * \frac{\partial U}{\partial t} + \frac{\partial E(U)}{\partial x} + \frac{\partial G(U)}{\partial y}  = S(U) \cr
 * U = \begin{bmatrix}\eta \cr q_x \cr q_y \end{bmatrix} \quad E = \begin{bmatrix}q_x \cr
 * \frac{q_x^2}{h} + \frac{1}{2}gh^2 \cr \frac{q_xq_y}{h}\end{bmatrix} \quad G = \begin{bmatrix}q_y \cr \frac{q_xq_y}{h} \cr
 * \frac{q_y^2}{h} + \frac{1}{2}gh^2 \end{bmatrix} \quad S = \begin{bmatrix}0 \cr -gh\frac{\partial z}{\partial x} \cr
 * -gh\frac{\partial z}{\partial y} \end{bmatrix}
 * \end{array} \end{equation} \f]
 *
 * Usages:
 * Use the 2 order basis with uniform mesh of 80 and 60 elements on x and y coordinate for test case ParabolicBowl :
 *
 *     mpirun -n 2 ./SWE2D 2 80 60 tri ParabolicBowl
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 * @todo
 *
 */

#include "SWEDriver2d.h"

void str2int(char *str, int *N, char* errmessage);

int main(int argc, char **argv){

    Solver *SWE2D = (Solver *)malloc(sizeof(Solver));

    int procid, nprocs; /* process number */

    MPI_Init(&argc, &argv);                 /* initialize MPI */
    MPI_Comm_rank(MPI_COMM_WORLD, &procid); /* read process id */
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); /* read process num */

    int Mx, My, N;
    /* get parameters */
    str2int(argv[1], &N, "Wrong degree input");
    str2int(argv[2], &Mx, "Wrong element number of x coordinate");
    str2int(argv[3], &My, "Wrong element number of y coordinate");
//    printf("N = %d, Mx = %d, My = %d\n", N, Mx, My);

    /* set stand element */
    StdRegions2d *shape;

    if ( !(memcmp(argv[4], "tri", 3)) ){
        shape = GenStdTriEle(N);
    }else if( !(memcmp(argv[4], "quad", 4)) ){
        shape = GenStdQuadEle(N);
    }else{
        printf("Wrong mesh type: %s\n", argv[4]);
        MPI_Finalize(); exit(1);
    }

    /* get mesh and physdomain */
    MultiReg2d   *mesh = SWEMesh2d(argv[5], shape, Mx, My);
    PhysDomain2d *phys = SWEInit2d(argv[5], mesh);


    FreeStdRegions2d(shape);
    FreeMultiReg2d(mesh);

    MPI_Finalize();
    return 0;
}

void str2int(char *str, int *N, char* errmessage){
    int info = sscanf(str,"%d",N);
    if (info!=1) {
        fprintf(stderr, "%s:%s \n", errmessage, str);
        exit(-1);
    }
}