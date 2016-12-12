#include "SWEDriver2d.h"
#include "LibUtilities/GenUniformMesh.h"

/* private variables */
MultiReg2d* ParabolicBowlMesh2d(char **argv, SWE_Solver2d *solver);
MultiReg2d* DamBreakMesh2d(char **argv, SWE_Solver2d *solver);
void GetCellNum_SWE2d(char **argv, int *Mx, int *My);
void GetOrder_SWE2d(char **argv, int *N);
StdRegions2d* AllocateStandShape_SWE2d(char **argv);
MultiReg2d* AllocateUniformMesh_SWE2d(char **argv, double xmin, double xmax, double ymin, double ymax);

/**
 * @brief
 * Setup mesh object and topography for various tests.
 *
 * @details
 * Setup the mesh object for different tests. The bottom topography
 * is stored in `solver` struct.
 * For test cases of `ParabolicBowl`, `DamBreakDry`, `DamBreakDry` and `FlowOver3Bumps`.
 * the mesh is uniform, and specified with # of cells on each coordinate.
 * For test cases of `TsunamiRunup`, the mesh is specified with mesh file name.
 *
 * @param[in]    argv   input argument.
 * @param[inout] solver SWE solver structure.
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * mesh     | MultiReg2d* | mesh object
 *
 * @todo
 * The user can specify either the mesh file name or the # of cell for uniform mesh.
 */
MultiReg2d* SWE_Mesh2d(char **argv, SWE_Solver2d *solver){

    MultiReg2d *mesh;
    /* setup mesh for various tests */
    if      ( !(memcmp(argv[1], "ParabolicBowl", 13)) ){
        mesh = ParabolicBowlMesh2d(argv, solver);
    }else if( !(memcmp(argv[1], "DamBreakDry"  , 11)) ){
        mesh = DamBreakMesh2d(argv, solver);
    }else if( !(memcmp(argv[1], "DamBreakWet"  , 11)) ){
        mesh = DamBreakMesh2d(argv, solver);
    }else{
        printf("Wrong name of test case: %s\n", argv[1]);
        MPI_Finalize(); exit(1);
    }

    return mesh;
}

/* allocate mesh for ParabolicBowl test case */
MultiReg2d* ParabolicBowlMesh2d(char **argv, SWE_Solver2d *solver){

    double xmin = -4000, xmax = 4000;
    double ymin = -4000, ymax = 4000;
    solver->gra = 9.81;
    double alpha = 1.6e-7;
    int Mx,My;
    /* # of cells on each coordinate */
    GetCellNum_SWE2d(argv, &Mx, &My);
    /* allocate uniform mesh grid */
    MultiReg2d *mesh = AllocateUniformMesh_SWE2d(argv, xmin, xmax, ymin, ymax);
    StdRegions2d *shape = mesh->stdcell;

    /* set bottom topography */
    int k, i;
    double r2;
    solver->bot = Matrix_create(mesh->K, shape->Np);
    for (k=0; k<mesh->K; k++){
        for (i=0; i<shape->Np; i++){
            r2 = mesh->x[k][i]*mesh->x[k][i] + mesh->y[k][i]*mesh->y[k][i];
            solver->bot[k][i] = alpha*r2;
        }
    }

    /* the grid length */
    solver->dx = min( (xmax - xmin)/(Mx+1.0)/shape->Nfp, (ymax - ymin)/(My+1)/shape->Nfp );
    return mesh;
}

/* allocate mesh for DamBreak test case */
MultiReg2d* DamBreakMesh2d(char **argv, SWE_Solver2d *solver){

    int Mx,My;
    double xmin = 0,    xmax = 1000;
    double ymin = -100, ymax = 100;
    solver->gra = 9.81;

    /* # of cells on each coordinate */
    GetCellNum_SWE2d(argv, &Mx, &My);
    /* allocate uniform mesh grid */
    MultiReg2d *mesh = AllocateUniformMesh_SWE2d(argv, xmin, xmax, ymin, ymax);
    StdRegions2d *shape = mesh->stdcell;

    /* set bottom topography */
    int k, i;
    solver->bot = Matrix_create(mesh->K, shape->Np);
    for (k=0; k<mesh->K; k++){
        for (i=0; i<shape->Np; i++){
            solver->bot[k][i] = 0.0;
        }
    }
    /* the grid length */
    solver->dx = min( (xmax - xmin)/(Mx+1.0)/shape->Nfp, (ymax - ymin)/(My+1)/shape->Nfp );
    return mesh;
}

/**
 * @brief
 * Allocate a uniform mesh from input arguments.
 *
 * @param[in] argv input arguments
 * @param[in] xmin minimal x coordinate
 * @param[in] xmax maximal x coordinate
 * @param[in] ymin minimal y coordinate
 * @param[in] ymax maximal y coordinate
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * mesh  | MultiReg2d* | pointer to the mesh structure
 *
 */
MultiReg2d* AllocateUniformMesh_SWE2d(char **argv, double xmin, double xmax, double ymin, double ymax){
    MultiReg2d *mesh;
    int Mx,My;

    GetCellNum_SWE2d(argv, &Mx, &My);
    StdRegions2d *shape = AllocateStandShape_SWE2d(argv);

    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    UnstructMesh *grid;
    if ( !(memcmp(argv[3], "tri", 3)) ){
        grid = ParallelUniformTriMesh_create(
                Mx, My, xmin, xmax, ymin, ymax, 1, procid, nprocs);
//        ParallelUniformTriMesh_create(Mx, My, xmin, xmax, ymin, ymax, 1,
//                                  procid, nprocs, &Ne, &Nv, &parEToV, &VX, &VY);
        mesh = MultiReg2d_create(shape, grid);
    }else if( !(memcmp(argv[3], "quad", 4)) ){
        grid = ParallelUniformQuadMesh_create(
                Mx, My, xmin, xmax, ymin, ymax, procid, nprocs);
//        ParallelUniformQuadMesh_create(Mx, My, xmin, xmax, ymin, ymax,
//                                   procid, nprocs, &Ne, &Nv, &parEToV, &VX, &VY);
        mesh = MultiReg2d_create(shape, grid);
    }else{
        printf("Wrong mesh type: %s.\n"
                       "The input should be either \'tri\' or \'quad\'.\n", argv[4]);
        MPI_Finalize(); exit(1);
    }

    UnstructMesh_free(grid);
    return mesh;
}

/* Allocate a standard element */
StdRegions2d* AllocateStandShape_SWE2d(char **argv){

    int N;
    GetOrder_SWE2d(argv, &N);

    StdRegions2d *shape;
    if ( !(memcmp(argv[3], "tri", 3)) ){
        shape = StdRegions2d_create(N, TRIANGLE);
    }else if( !(memcmp(argv[3], "quad", 4)) ){
        shape = StdRegions2d_create(N, QUADRIL);
    }else{
        printf("Wrong mesh type: %s.\n"
                       "The input should be either \'tri\' or \'quad\'.\n", argv[3]);
        MPI_Finalize(); exit(1);
    }
    return shape;
}

/* get the input degree */
void GetOrder_SWE2d(char **argv, int *N){
    str2int(argv[2], N, "Wrong degrees");
    return;
}

/* the # of cell for uniform mesh grid */
void GetCellNum_SWE2d(char **argv, int *Mx, int *My){
    str2int(argv[4], Mx, "Wrong element number of x coordinate");
    str2int(argv[5], My, "Wrong element number of y coordinate");
    return;
}