#include "SWEDriver2d.h"

MultiReg2d* ParabolicBowlMesh2d(SWE_Solver2d *solver, StdRegions2d *shape, int Mx, int My);
MultiReg2d* DamBreakMesh2d     (SWE_Solver2d *solver, StdRegions2d *shape, int Mx, int My);

/**
 * @brief
 * Setup mesh domain and topography for various tests.
 *
 * @details
 * Setup the mesh domain for different tests. The number of elements is defined by `Mx` and `My`.
 * The bottom topography is stored in `solver` struct.
 *
 * @param[in] casename  name of test case
 * @param[in] shape     standard element
 * @param[in] Mx        number of elements on x coordinate
 * @param[in] My        number of elements on y coordinate
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * mesh     | MultiReg2d* | mesh object
 *
 */
MultiReg2d* SWE_Mesh2d(char *casename, SWE_Solver2d *solver, StdRegions2d *shape, int Mx, int My){
    MultiReg2d *mesh;

    /* setup mesh for various tests */
    if      ( !(memcmp(casename, "ParabolicBowl", 13)) ){
        mesh = ParabolicBowlMesh2d(solver, shape, Mx, My);
    }else if( !(memcmp(casename, "DamBreakDry"  , 11)) ){
        mesh = DamBreakMesh2d(solver, shape, Mx, My);
    }else if( !(memcmp(casename, "DamBreakWet"  , 11)) ){
        mesh = DamBreakMesh2d(solver, shape, Mx, My);
    }else{
        printf("Wrong name of test case: %s\n", casename);
        MPI_Finalize(); exit(1);
    }

    return mesh;
}

MultiReg2d* ParabolicBowlMesh2d(SWE_Solver2d *solver, StdRegions2d *shape, int Mx, int My){
    MultiReg2d *mesh;

    double xmin = -4000, xmax = 4000;
    double ymin = -4000, ymax = 4000;
    solver->gra = 9.81;
    double alpha = 1.6e-7;
    double w    = sqrt(8*solver->gra*alpha);

    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int Ne, Nv;
    int **parEToV;
    double *VX, *VY;

    /* generate triangle or quadrilateral mesh */
    switch (shape->Nv) {
        case 3:
            GenParallelUniformTriMesh(Mx, My, xmin , xmax, ymin, ymax, 1,
                                      procid, nprocs, &Ne, &Nv, &parEToV, &VX, &VY);
            mesh = GenMultiReg2d(shape, Ne, Nv, parEToV, VX, VY);
            break;
        case 4:
            GenParallelUniformQuadMesh(Mx, My, xmin, xmax, ymin, ymax,
                                       procid, nprocs, &Ne, &Nv, &parEToV, &VX, &VY);
            mesh = GenMultiReg2d(shape, Ne, Nv, parEToV, VX, VY);
            break;
        default:
            printf("Mesh Generator error! Wrong number of element vertex: %d\n", shape->Nv);
            MPI_Finalize(); exit(2);
    }

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

    DestroyIntMatrix(parEToV);
    DestroyVector(VX);
    DestroyVector(VY);

    /* the grid length */
    solver->dx = min( (xmax - xmin)/(Mx+1.0)/shape->Nfp, (ymax - ymin)/(My+1)/shape->Nfp );

    return mesh;
}

MultiReg2d* DamBreakMesh2d(SWE_Solver2d *solver, StdRegions2d *shape, int Mx, int My){
    MultiReg2d *mesh;

    double xmin = 0,    xmax = 1000;
    double ymin = -100, ymax = 100;
    solver->gra = 9.81;

    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int Ne, Nv;
    int **parEToV;
    double *VX, *VY;

    /* generate triangle or quadrilateral mesh */
    switch (shape->Nv) {
        case 3:
            GenParallelUniformTriMesh(Mx, My, xmin , xmax, ymin, ymax, 1,
                                      procid, nprocs, &Ne, &Nv, &parEToV, &VX, &VY);
            mesh = GenMultiReg2d(shape, Ne, Nv, parEToV, VX, VY);
            break;
        case 4:
            GenParallelUniformQuadMesh(Mx, My, xmin, xmax, ymin, ymax,
                                       procid, nprocs, &Ne, &Nv, &parEToV, &VX, &VY);
            mesh = GenMultiReg2d(shape, Ne, Nv, parEToV, VX, VY);
            break;
        default:
            printf("Mesh Generator error! Wrong number of element vertex: %d\n", shape->Nv);
            MPI_Finalize(); exit(2);
    }

    /* set bottom topography */
    int k, i;
    solver->bot = BuildMatrix(mesh->K, shape->Np);
    for (k=0; k<mesh->K; k++){
        for (i=0; i<shape->Np; i++){
            solver->bot[k][i] = 0.0;
        }
    }

    DestroyIntMatrix(parEToV);
    DestroyVector(VX);
    DestroyVector(VY);

    /* the grid length */
    solver->dx = min( (xmax - xmin)/(Mx+1.0)/shape->Nfp, (ymax - ymin)/(My+1)/shape->Nfp );

    return mesh;
}