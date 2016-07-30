#include "MultiRegionsTest.h"

int main(int argc, char **argv){
    int info;

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    int TriCoordTest(void);
    info = TriCoordTest();

    int QuadCoordTest(void);
    info = QuadCoordTest();

    MPI_Finalize();
    return info;
}

int QuadCoordTest(void){
    int i,j, N=3, K=2;
    int Np = (N+1)*(N+1), info;

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

    StdRegions2d *quad = GenStdQuadEle(N);
    MultiReg2d *mesh;
    SetTestQuadMesh(quad, mesh);

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
    int i,j, N=3;
    int K=4, Np = (N+1)*(N+2)/2;
    int info;

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

    StdRegions2d *tri = GenStdTriEle(N);
    MultiReg2d *mesh;
    SetTestTriMesh(tri, mesh);

    for(i=0;i<K;i++){
        for(j=0;j<Np;j++){
            temp[i*Np+j] = mesh->x[i][j];
            temp_exact[i*Np+j] = xt[j][i];
        }
    }
    info = CreateVectorTest("TriCoor x", temp, temp_exact, K*Np);

    for(i=0;i<K;i++){
        for(j=0;j<Np;j++){
            temp[i*Np+j] = mesh->y[i][j];
            temp_exact[i*Np+j] = yt[j][i];
        }
    }
    info = CreateVectorTest("TriCoor y", temp, temp_exact, K*Np);

    return info;
}