#include "SWEDriver2d.h"

MultiReg2d* ParabolicBowlMesh2d(char **argv, SWE_Solver2d *solver);
MultiReg2d* DamBreakMesh2d(char **argv, SWE_Solver2d *solver);
void        GetCellNum(char **argv, int *Mx, int *My);
void        GetOrder(char **argv, int *N);
StdRegions2d* AllocateStandShape(char **argv);
MultiReg2d*   AllocateUniformMesh(char **argv, double xmin, double xmax, double ymin, double ymax);

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
 * @param[char**]        argv   input argument.
 * @param[SWE_Solver2d*] solver SWE solver structure.
 *
 * @return
 * return values:
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
    GetCellNum(argv, &Mx, &My);
    /*  */
    MultiReg2d *mesh = AllocateUniformMesh(argv,xmin,xmax,ymin,ymax);
    StdRegions2d *shape = mesh->stdcell;

    /* set bottom topography */
    int k, i;
    double r2;
    solver->bot = BuildMatrix(mesh->K, shape->Np);
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
    GetCellNum(argv, &Mx, &My);
    MultiReg2d *mesh = AllocateUniformMesh(argv,xmin,xmax,ymin,ymax);
    StdRegions2d *shape = mesh->stdcell;

    /* set bottom topography */
    int k, i;
    solver->bot = BuildMatrix(mesh->K, shape->Np);
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
 * @param[char**] argv input arguments
 * @param[double] xmin minimal x coordinate
 * @param[double] xmax maximal x coordinate
 * @param[double] ymin minimal y coordinate
 * @param[double] ymax maximal y coordinate
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * mesh  | MultiReg2d* | pointer to the mesh structure
 *
 */
MultiReg2d* AllocateUniformMesh(char **argv, double xmin, double xmax, double ymin, double ymax){
    MultiReg2d *mesh;
    int Ne,Nv,Mx,My;
    int **parEToV;
    double *VX,*VY;

    GetCellNum(argv, &Mx, &My);
    StdRegions2d *shape = AllocateStandShape(argv);

    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if ( !(memcmp(argv[3], "tri", 3)) ){
        GenParallelUniformTriMesh(Mx, My, xmin, xmax, ymin, ymax, 1,
                                  procid, nprocs, &Ne, &Nv, &parEToV, &VX, &VY);
        mesh = GenMultiReg2d(shape, Ne, Nv, parEToV, VX, VY);
    }else if( !(memcmp(argv[3], "quad", 4)) ){
        GenParallelUniformQuadMesh(Mx, My, xmin, xmax, ymin, ymax,
                                   procid, nprocs, &Ne, &Nv, &parEToV, &VX, &VY);
        mesh = GenMultiReg2d(shape, Ne, Nv, parEToV, VX, VY);
    }else{
        printf("Wrong mesh type: %s.\n"
                       "The input should be either \'tri\' or \'quad\'.\n", argv[4]);
        MPI_Finalize(); exit(1);
    }

    DestroyIntMatrix(parEToV);
    DestroyVector(VX);
    DestroyVector(VY);
    return mesh;
}

/**
 * @brief
 * Allocate the standard element.
 *
 * @param[char**] argv input arguments.
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * shape   | StdRegions2d* | the standard element structure
 *
 */
StdRegions2d* AllocateStandShape(char **argv){
    StdRegions2d *shape;
    int N;
    GetOrder(argv, &N);

    if ( !(memcmp(argv[3], "tri", 3)) ){
        shape = GenStdTriEle(N);
    }else if( !(memcmp(argv[3], "quad", 4)) ){
        shape = GenStdQuadEle(N);
    }else{
        printf("Wrong mesh type: %s.\n"
                       "The input should be either \'tri\' or \'quad\'.\n", argv[4]);
        MPI_Finalize(); exit(1);
    }
    return shape;
}

/* get the input degree */
void GetOrder(char **argv, int *N){
    str2int(argv[2], N, "Wrong degrees");
    return;
}

/* the # of cell for uniform mesh grid */
void GetCellNum(char **argv, int *Mx, int *My){
    str2int(argv[4], Mx, "Wrong element number of x coordinate");
    str2int(argv[5], My, "Wrong element number of y coordinate");
    return;
}