#include "MultiRegionsTest.h"

int main(int argc, char **argv){

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    int info;

    int TriCoordTest(void);
    info = TriCoordTest();

    int QuadCoordTest(void);
    info = QuadCoordTest();

    MPI_Finalize();
    return info;
}

int QuadCoordTest(void){
    int i,j, N=3, Nvert=4;
    int K=2, Nv = 6;
    int Np = (N+1)*(N+1);
    int info;

    int **EToV = BuildIntMatrix(K, Nvert);
    double *VX = BuildVector(Nv);
    double *VY = BuildVector(Nv);

    int EToVt[2][4] = {{1,4,5,2},{2,5,6,3}};
    double VXt[6] = {-1, 0, 1, -1, 0, 1};
    double VYt[6] = {1, 1, 1, 0, 0, 0};

    double xt[16][2] = {{-1.000000e+00,0.000000e+00,},
                        {-1.000000e+00,0.000000e+00,},
                        {-1.000000e+00,0.000000e+00,},
                        {-1.000000e+00,0.000000e+00,},
                        {-7.236068e-01,2.763932e-01,},
                        {-7.236068e-01,2.763932e-01,},
                        {-7.236068e-01,2.763932e-01,},
                        {-7.236068e-01,2.763932e-01,},
                        {-2.763932e-01,7.236068e-01,},
                        {-2.763932e-01,7.236068e-01,},
                        {-2.763932e-01,7.236068e-01,},
                        {-2.763932e-01,7.236068e-01,},
                        {0.000000e+00,1.000000e+00,},
                        {0.000000e+00,1.000000e+00,},
                        {0.000000e+00,1.000000e+00,},
                        {0.000000e+00,1.000000e+00,}};
    double yt[16][2] = {{1.000000e+00,1.000000e+00,},
                        {7.236068e-01,7.236068e-01,},
                        {2.763932e-01,2.763932e-01,},
                        {0.000000e+00,0.000000e+00,},
                        {1.000000e+00,1.000000e+00,},
                        {7.236068e-01,7.236068e-01,},
                        {2.763932e-01,2.763932e-01,},
                        {0.000000e+00,0.000000e+00,},
                        {1.000000e+00,1.000000e+00,},
                        {7.236068e-01,7.236068e-01,},
                        {2.763932e-01,2.763932e-01,},
                        {0.000000e+00,0.000000e+00,},
                        {1.000000e+00,1.000000e+00,},
                        {7.236068e-01,7.236068e-01,},
                        {2.763932e-01,2.763932e-01,},
                        {0.000000e+00,0.000000e+00,}};

    double **xt_ext = BuildMatrix(K, Np);
    double **yt_ext = BuildMatrix(K, Np);

    for(i=0;i<Nv;i++){
        VX[i] = VXt[i];
        VY[i] = VYt[i];
    }

    for(i=0;i<K;i++){
        for(j=0;j<Nvert;j++) {
            EToV[i][j] = EToVt[i][j] - 1;
        }
    }


    StdRegions2d *quad = GenStdQuadEle(N);
    MultiReg2d *mesh = GenMultiReg2d(quad, K, Nv, EToV, VX, VY);

    for(i=0;i<K;i++){
        for(j=0;j<quad->Np;j++){
            xt_ext[i][j] = xt[j][i];
        }
    }

    info = CreateMatrixTest("quad x", mesh->x, xt_ext, K, Np);

    for(i=0;i<K;i++){
        for(j=0;j<Np;j++){
            yt_ext[i][j] = yt[j][i];
        }
    }

    info = CreateMatrixTest("quad y", mesh->y, yt_ext, K, Np);

    return info;
}

int TriCoordTest(void){
    int i,j, N=3, Nvert=3;
    int K=4, Nv = 6;
    int Np = (N+1)*(N+2)/2;
    int info;

    int **EToV = BuildIntMatrix(K, Nvert);
    double *VX = BuildVector(Nv);
    double *VY = BuildVector(Nv);

    int EToVt[4][3] = {{1,4,2},{2,4,5},
                       {2,5,6},{2,6,3}};
    double VXt[6] = {-1, 0, 1, -1, 0, 1};
    double VYt[6] = {1, 1, 1, 0, 0, 0};

    double xt[10][4] = {{-1,	0,	0,	0},
                        {-1,	-0.276393202250021,	5.55111512312578e-17,	0.276393202250021},
                        {-1,	-0.723606797749979,	0,	0.723606797749979},
                        {-1,	-1,	0,	1},
                        {-0.723606797749979,	1.11022302462516e-16,	0.276393202250021,	0.276393202250021},
                        {-0.666666666666667,	-0.333333333333333,	0.333333333333333,	0.666666666666667},
                        {-0.723606797749979,	-0.723606797749979,	0.276393202250021,	1},
                        {-0.276393202250021,	0,	0.723606797749979,	0.723606797749979},
                        {-0.276393202250021,	-0.276393202250021,	0.723606797749979,	1},
                        {0,	0,	1,	1}};
    double yt[10][4] = {{1,	1,	1,	1},
                        {0.723606797749979,	0.723606797749979,	0.723606797749979,	0.723606797749979},
                        {0.276393202250021,	0.276393202250021,	0.276393202250021,	0.276393202250021},
                        {0,	0,	0,	0},
                        {1,	0.723606797749979,	0.723606797749979,	1},
                        {0.666666666666667,	0.333333333333333,	0.333333333333333,	0.666666666666667},
                        {0.276393202250021,	0,	0,	0.276393202250021},
                        {1,	0.276393202250021,	0.276393202250021,	1},
                        {0.723606797749979,	5.55111512312578e-17,	5.55111512312578e-17,	0.723606797749979},
                        {1,	0,	0,	1}};

    double temp[40], temp_exact[40];

    for(i=0;i<Nv;i++){
        VX[i] = VXt[i];
        VY[i] = VYt[i];
    }

    for(i=0;i<K;i++){
        for(j=0;j<Nvert;j++)
            EToV[i][j] = EToVt[i][j]-1;
    }

    StdRegions2d *tri = GenStdTriEle(N);
    MultiReg2d *mesh = GenMultiReg2d(tri, K, Nv, EToV, VX, VY);

    for(i=0;i<K;i++){
        for(j=0;j<Np;j++){
            temp[i*Np+j] = mesh->x[i][j];
            temp_exact[i*Np+j] = xt[j][i];
        }
    }

    info = CreateVectorTest("tri x", temp, temp_exact, K*Np);

    for(i=0;i<K;i++){
        for(j=0;j<Np;j++){
            temp[i*Np+j] = mesh->y[i][j];
            temp_exact[i*Np+j] = yt[j][i];
        }
    }

    info = CreateVectorTest("tri y", temp, temp_exact, K*Np);

    return info;
}