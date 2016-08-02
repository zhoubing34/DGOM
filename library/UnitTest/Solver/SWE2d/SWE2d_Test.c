//
// Created by li12242 on 16/8/2.
//

#include <PhysDomain/PhysDomain.h>
#include "UnitTest/MultiRegions/MultiRegionsTest.h"
#include "SWE2d/SWEDriver2d.h"

int main(int argc, char **argv){

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    /* set mesh and physical domain */
    int N=1, Nfields=3;
    StdRegions2d *tri = GenStdTriEle(N);
    MultiReg2d   *triMesh;
    SetTestTriMesh(tri, triMesh);
    PhysDomain2d *phys = GenPhysDomain2d(triMesh, Nfields);
    int i,k;
    for(k=0;k<triMesh->K;k++){
        for(i=0;i<tri->Np;i++) {
            int ind = (k*tri->Np + i)*Nfields;
            /* initial value assignment */
            if(k<=2){
                phys->f_Q[ind++] = 2.0;
            }else{
                phys->f_Q[ind++] = 0.0;
            }
            phys->f_Q[ind++] = 0;
            phys->f_Q[ind++] = 0;
        }
    }

    /* solver */
    SWESolver *solver = (SWESolver*) calloc(1, sizeof(SWESolver));
    solver->hcrit = 1.0e-4;
    solver->bot   = BuildMatrix(triMesh->K, tri->Np); /* set bottom topography */

    for(k=0;k<triMesh->K;k++){
        for(i=0;i<tri->Np;i++)
            solver->bot[k][i] = triMesh->x[k][i];
    }

    return 0;
}

void FluxTest(){}