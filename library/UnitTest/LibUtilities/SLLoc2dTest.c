#include <LibUtilities/SlopeLimiter.h>
#include <UnitTest/MultiRegions/MultiRegionsTest.h>

void PrintVec(char *message, real *vec, int n);
void PrintVar(PhysDomain2d *phys);

int main(int argc, char **argv){

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    /* set mesh and physical domain */
    int N=1, Nfields=3;
    StdRegions2d *tri = StdTriEle_create(N);
    MultiReg2d   *mesh;
    SetTestTriMesh(tri, mesh);
    PhysDomain2d *phys = GenPhysDomain2d(mesh, Nfields);
    int i,k;
    for(k=0;k<mesh->K;k++){
        for(i=0;i<tri->Np;i++) {
            int ind = (k*tri->Np + i)*Nfields;
            /* initial value assignment */
            if(k<=1){
                phys->f_Q[ind++] = (real)mesh->x[k][i]+1;
                phys->f_Q[ind++] = (real)mesh->x[k][i]+1;
                phys->f_Q[ind++] = (real)mesh->x[k][i]+1;
            }else{
                phys->f_Q[ind++] = 0.0;
                phys->f_Q[ind++] = 0.0;
                phys->f_Q[ind++] = 0.0;
            }
        }
    }

    for(i=0;i<5;i++) {
        if(mesh->procid==0) {
            printf("Initial conditions\n");
            PrintVar(phys);
        }
        SLLoc2d(phys, 1.0);

        if(mesh->procid==0) {
            printf("Initial conditions\n");
            PrintVar(phys);
        }
    }

    MPI_Finalize();
    return 0;
}

void PrintVar(PhysDomain2d *phys){
    MultiReg2d   *mesh  = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    real *var = (real *) malloc(mesh->K*shape->Np*sizeof(real));
    int k,i,ind,sk=0;
    int Np = shape->Np, Nfields=phys->Nfields;
    for (k=0;k<mesh->K;k++){
        for (i=0;i<Np;i++){
            ind = (k*Np + i)*Nfields;
            // printf("k = %d, i = %d, ind = %d\n", k, i, ind);
            var[sk++] = phys->f_Q[ind];
        }
        PrintVec("h ",var+k*Np, Np);
    }

    sk = 0;
    for (k=0;k<mesh->K;k++){
        for (i=0;i<Np;i++){
            ind = (k*Np + i)*Nfields + 1;
            // printf("k = %d, i = %d, ind = %d\n", k, i, ind);
            var[sk++] = phys->f_Q[ind];
        }
        PrintVec("qx",var+k*Np, Np);
    }

    sk = 0;
    for (k=0;k<mesh->K;k++){
        for (i=0;i<Np;i++){
            ind = (k*Np + i)*Nfields + 2;
            // printf("k = %d, i = %d, ind = %d\n", k, i, ind);
            var[sk++] = phys->f_Q[ind];
        }
        PrintVec("qy",var+k*Np, Np);
    }

    free(var);
}

void PrintVec(char *message, real *vec, int n){
    int i;
    printf("%s: ", message);
    for (i=0;i<n;i++){
        printf("%f, ", vec[i]);
    }
    printf("\n");
}
