//
// Created by li12242 on 16/6/12.
//

#include "test.h"

#define N 5
#define Np (N+1)*(N+2)/2
int main(int argc, char **argv){

    void GetTriCoordTest(void);
    GetTriCoordTest();
    return 0;
}


void GetTriCoordTest(void){
    double *r = BuildVector(Np);
    double *s = BuildVector(Np);
    int i;

    GetTriCoord(N, r, s);

    for(i=0;i<Np;i++){
        printf("r[%d]=%12.4f, s[%d]=%12.4f\n", i,r[i],i,s[i]);
    }

    DestroyVector(r);
    DestroyVector(s);

}

void WarpfactorTest(void){
    double *v = BuildVector(Np);
    double *w = BuildVector(Np);
    int i;

    for(i=0;i<Np;i++){
        v[i] = i/(Np-1.0)*2.0 - 1.0;
        printf("v[%d] = %lg,\t w[%d] = %lg\n", i, v[i], i, w[i]);
    }

    Warpfactor(N, v, Np, w);

    for(i=0;i<Np;i++){
        printf("v[%d] = %lg,\t w[%d] = %e\n", i, v[i], i, w[i]);
    }

    DestroyVector(v);
}
