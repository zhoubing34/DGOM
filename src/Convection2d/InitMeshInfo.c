#include "Convection2d/Convection2d.h"

/* get the surface  */
double InitTriMeshInfo(Mesh * mesh, int Nfields){
    printf("Np = %d, BSIZE = %d\n", p_Np, BSIZE);


    /*  float LIFT  */
    int sz = p_Np*(p_Nfp)*(p_Nfaces)*sizeof(float);
    mesh->f_LIFT = (float*) malloc(sz);
    int sk = 0, n, m, f, k;

    for(n=0;n<p_Np;++n){
        for(m=0;m<p_Nfp*p_Nfaces;++m){
            mesh->f_LIFT[sk++] = (float) mesh->LIFT[n][m];
        }
    }

    /*  float Dr & Ds */
    sz = p_Np*p_Np*sizeof(float);
    mesh->f_Dr = (float*) malloc(sz);
    mesh->f_Ds = (float*) malloc(sz);

    sk = 0;
    for(n=0;n<p_Np;++n){
        for(m=0;m<p_Np;++m){
            mesh->f_Dr[sk] = mesh->Dr[n][m];
            mesh->f_Ds[sk] = mesh->Ds[n][m];
            ++sk;
        }
    }

    /* vgeo */
    double drdx, dsdx, drdy, dsdy, J;
    mesh->vgeo = (float*) calloc(4*mesh->K, sizeof(float));

    for(k=0;k<mesh->K;++k){
        GeometricFactorsTri(mesh, k, &drdx, &dsdx, &drdy, &dsdy, &J);
        mesh->vgeo[k*4+0] = drdx;
        mesh->vgeo[k*4+1] = drdy;
        mesh->vgeo[k*4+2] = dsdx;
        mesh->vgeo[k*4+3] = dsdy;
    }

    /* surfinfo (vmapM, vmapP, Fscale, Bscale, nx, ny, nz, 0) */
    sz = mesh->K*p_Nfp*p_Nfaces*6*sizeof(float);
    mesh->surfinfo = (float*) malloc(sz);

    /* local-local info */
    sk = 0;
    int skP = -1;
    double *nxk = BuildVector(mesh->Nfaces);
    double *nyk = BuildVector(mesh->Nfaces);
    double *sJk = BuildVector(mesh->Nfaces);

    double dt = 1e6;

    sk = 0;
    for(k=0;k<mesh->K;++k){
        GeometricFactorsTri(mesh, k, &drdx, &dsdx, &drdy, &dsdy, &J);
        NormalsTri(mesh, k, nxk, nyk, sJk);

        for(f=0;f<mesh->Nfaces;++f){
            dt = min(dt, J/sJk[f]);

            for(m=0;m<p_Nfp;++m){

                int id  = m + f*p_Nfp + p_Nfp*p_Nfaces*k;
                int idM = mesh->vmapM[id];
                int idP = mesh->vmapP[id];
                int  nM = idM%p_Np;
                int  nP = idP%p_Np;
                int  kM = (idM-nM)/p_Np;
                int  kP = (idP-nP)/p_Np;

                idM = Nfields*(nM+p_Np*kM);
                idP = Nfields*(nP+p_Np*kP);

                /* stub resolve some other way */
                if(mesh->vmapP[id]<0){
                    idP = mesh->vmapP[id]; /* -ve numbers */
                }

                mesh->surfinfo[sk++] = idM;
                mesh->surfinfo[sk++] = idP;
                mesh->surfinfo[sk++] = sJk[f]/(2*J);
                mesh->surfinfo[sk++] = (idM==idP)?-1.:1.;
                mesh->surfinfo[sk++] = nxk[f];
                mesh->surfinfo[sk++] = nyk[f];
            }
        }
    }

    return dt;
}