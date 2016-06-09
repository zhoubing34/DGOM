#include "Convection2d/Convection2d.h"

/**
 * @brief
 *
 *
 * @details
 *
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 * @return
 * return values：
 * name     | type     | description of value
 * -------- |----------|----------------------
 * car_id   | int      |
 * car_info | object   |
 *
 */

double InitMeshInfo(Mesh * mesh, int Nfields){

    int procid = mesh->procid;
    int K = mesh->K;

#if defined DEBUG
    printf("Procs %d: Entering InitMeshInfo\n", procid);
#endif

    printf("Np = %d, BSIZE = %d\n", p_Np, BSIZE);


    /* allocate memory for variables */
    mesh->f_Q    = (float*) calloc(K*BSIZE*p_Nfields, sizeof(float));
    mesh->f_rhsQ = (float*) calloc(K*BSIZE*p_Nfields, sizeof(float));
    mesh->f_resQ = (float*) calloc(K*BSIZE*p_Nfields, sizeof(float));

    mesh->f_s = (float*) calloc(K*BSIZE*2, sizeof(float));

    mesh->tcflag    = (float *) calloc(K, sizeof(float));
    mesh->area      = (float *) calloc(K, sizeof(float));
    mesh->ciradius  = (float*) calloc(K, sizeof(float));

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

    /* volume geometric factor and
     * radius circumscribed circle  */
    double *drdx, *dsdx, *drdy, *dsdy, *J;
    sz = p_Np*sizeof(double);

    drdx = (double*) malloc(sz);
    drdy = (double*) malloc(sz);
    dsdx = (double*) malloc(sz);
    dsdy = (double*) malloc(sz);
    J    = (double*) malloc(sz);

    mesh->vgeo = (float*) calloc(4*mesh->K*p_Np, sizeof(float));
    mesh->J    = (float*) calloc(K*p_Np, sizeof(float));
    int sj = 0;
    sk = 0;
    for(k=0;k<mesh->K;++k){
        GeometricFactors(mesh, k, drdx, dsdx, drdy, dsdy, J);
        for(n=0;n<p_Np;n++){
            mesh->vgeo[sk++] = (float) drdx[n];
            mesh->vgeo[sk++] = (float) drdy[n];
            mesh->vgeo[sk++] = (float) dsdx[n];
            mesh->vgeo[sk++] = (float) dsdy[n];
//            mesh->vgeo[sk++] = J[n];
            mesh->J[sj++] = (float) J[n];

            mesh->area[k] += (float) J[n]*mesh->wv[n];
        }
        mesh->ciradius[k] = sqrt(mesh->area[k]/M_PI);
    }

    /* surfinfo (vmapM, vmapP, Fscale, Bscale, nx, ny, nz, 0) */
    sz = mesh->K*p_Nfp*p_Nfaces*6*sizeof(float);
    mesh->surfinfo = (float*) malloc(sz);

    /* local-local info */
    double *nxk = BuildVector(mesh->Nfaces);
    double *nyk = BuildVector(mesh->Nfaces);
    double *sJk = BuildVector(mesh->Nfaces);

    double dt = 1e6;

    sk = 0;
    for(k=0;k<mesh->K;++k){
        GeometricFactors(mesh, k, drdx, dsdx, drdy, dsdy, J);
        Normals(mesh, k, nxk, nyk, sJk);

        for(f=0;f<mesh->Nfaces;++f){

            for(m=0;m<p_Nfp;++m){
                dt = min(dt, J[mesh->Fmask[f][m]]/sJk[f]);

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
                mesh->surfinfo[sk++] = sJk[f]/(2*J[mesh->Fmask[f][m]]);
                mesh->surfinfo[sk++] = sJk[f];
                mesh->surfinfo[sk++] = nxk[f];
                mesh->surfinfo[sk++] = nyk[f];
            }
        }
    }

    /* deallocate mem */
    free(drdx);
    free(drdy);
    free(dsdx);
    free(dsdy);
    free(J);

    DestroyVector(nxk);
    DestroyVector(nyk);
    DestroyVector(sJk);

#if defined DEBUG
    printf("Procs %d: Leaving InitMeshInfo\n", procid);
#endif

    return dt;
}