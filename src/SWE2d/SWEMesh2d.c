#include "SWEDriver2d.h"

MultiReg2d* ParabolicBowlInit(StdRegions2d *shape, int Mx, int My);
MultiReg2d* DamBreakDry      (StdRegions2d *shape, int Mx, int My);

/**
 * @brief
 * Setup mesh domain for various tests.
 *
 * @details
 * Setup the mesh domain for different tests. The number of elements is defined by `Mx` and `My`.
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
 *
 */
MultiReg2d* SWEMesh2d(char *casename, StdRegions2d *shape, int Mx, int My){
    MultiReg2d *mesh;

    /* setup mesh for various tests */
    if      ( !(memcmp(casename, "ParabolicBowl", 13)) ){
        mesh = ParabolicBowlInit(shape, Mx, My);
    }else if( !(memcmp(casename, "DamBreakDry"  , 11)) ){
        mesh = DamBreakDry(shape, Mx, My);
    }else{
        printf("Wrong name of test case: %s\n", casename);
        MPI_Finalize(); exit(1);
    }

    return mesh;
}

MultiReg2d* ParabolicBowlInit(StdRegions2d *shape, int Mx, int My){
    MultiReg2d *mesh;

    return mesh;
}

MultiReg2d* DamBreakDry(StdRegions2d *shape, int Mx, int My){
    MultiReg2d *mesh;

    return mesh;
}