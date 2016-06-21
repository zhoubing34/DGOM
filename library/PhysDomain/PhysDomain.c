/**
 * @file
 * Set physical variables
 *
 * @brief
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "PhysDomain.h"

/* public variables */
void SetSurfInfo2d(MultiReg2d *mesh, int Nfields, float *surfinfo);
void SetMapOut2d(MultiReg2d *mesh, int Nfields, int *parmapOUT);


PhysDomain2d* GetPhysDomain2d(MultiReg2d *mesh, int Nfields){

    PhysDomain2d *phys = (PhysDomain2d *) calloc(1, sizeof(PhysDomain2d));
    StdRegions2d *shape = mesh->stdcell;

    int K  = mesh->K;   /* number of elements */
    int Np = shape->Np; /* number of nodes in each element */

    int Nfp = shape->Nfp;
    int Nfaces = shape->Nfaces;

    /* number of variables */
    phys->Nfields = Nfields;

    /* mesh */
    phys->mesh = mesh;

    /* volume info */
    phys->vgeo = mesh->vgeo;

    /* data array */
    phys->f_Q    = (float *) calloc(K*Np*Nfields, sizeof(float));
    phys->f_rhsQ = (float *) calloc(K*Np*Nfields, sizeof(float));
    phys->f_resQ = (float *) calloc(K*Np*Nfields, sizeof(float));

    /* MPI send/recv buffer */
    phys->f_inQ  = (float *) calloc(mesh->parNtotalout*Nfields, sizeof(float));
    phys->f_outQ = (float *) calloc(mesh->parNtotalout*Nfields, sizeof(float));

    phys->parNtotalout = mesh->parNtotalout*Nfields;
    phys->parmapOUT = BuildIntVector(phys->parNtotalout);
    SetMapOut2d(mesh, Nfields, phys->parmapOUT);

    /* surface info */
    int sz = K*Nfp*Nfaces*6*sizeof(float);
    phys->surfinfo = (float*) malloc(sz);
    SetSurfInfo2d(mesh, Nfields, phys->surfinfo);

    return phys;
}

void FreePhysDomain2d(PhysDomain2d *phys){
    free(phys->f_Q);
    free(phys->f_resQ);
    free(phys->f_rhsQ);
    free(phys->f_inQ);
    free(phys->f_outQ);

    DestroyIntVector(phys->parmapOUT);
    free(phys->surfinfo);
}

void SetMapOut2d(MultiReg2d *mesh, int Nfields, int *parmapOUT){
    int p2, n1, m, fld;
    int nprocs = mesh->nprocs;
    int procid = mesh->procid;
    StdRegions2d *shape = mesh->stdcell;
    int Nfp = shape->Nfp;
    int Np  = shape->Np;

    int sp=0, sk=0;
    for(p2=0;p2<nprocs;++p2){
        if(p2!=procid) {
            /* for each received face */
            for (m = 0; m < mesh->Npar[p2]; ++m) {
                for (n1 = 0; n1 < Nfp; ++n1) {
                    for (fld = 0; fld < Nfields; ++fld) {
                        /* sk and node index determine the map relationship */
                        parmapOUT[sp++] = Nfields * (mesh->parmapOUT[sk]) + fld;
                    }
                    sk++;
                }
            }
        }
    }
}

void SetSurfInfo2d(MultiReg2d *mesh, int Nfields, float *surfinfo){

    int k,f,m,sk = 0;
    StdRegions2d *shape = mesh->stdcell;

    int K  = mesh->K;
    int Np = shape->Np;
    int Nfaces = shape->Nfaces;
    int Nfp = shape->Nfp;

    double *drdx, *dsdx, *drdy, *dsdy, *J;
    int sz = Np*sizeof(double);

    drdx = (double*) malloc(sz);
    drdy = (double*) malloc(sz);
    dsdx = (double*) malloc(sz);
    dsdy = (double*) malloc(sz);
    J    = (double*) malloc(sz);

    /* local-local info */
    double *nxk = BuildVector(Nfaces);
    double *nyk = BuildVector(Nfaces);
    double *sJk = BuildVector(Nfaces);

    for(k=0;k<mesh->K;++k){
        GeoFactor2d(shape->Np, mesh->x[k], mesh->y[k],
                    shape->Dr, shape->Ds,
                    drdx, dsdx, drdy, dsdy, J);
        Normals2d(Nfaces, mesh->GX[k], mesh->GY[k], nxk, nyk, sJk);

        for(f=0;f<Nfaces;++f){

            for(m=0;m<Nfp;++m){

                int id  = m + f*Nfp + Nfp*Nfaces*k;
                int idM = mesh->vmapM[id];
                int idP = mesh->vmapP[id];
                int  nM = idM%Np;
                int  nP = idP%Np;
                int  kM = (idM-nM)/Np;
                int  kP = (idP-nP)/Np;

                idM = Nfields*(nM+Np*kM);
                idP = Nfields*(nP+Np*kP);

                /* stub resolve some other way */
                if(mesh->vmapP[id]<0){
                    idP = mesh->vmapP[id]; /* -ve numbers */
                }

                surfinfo[sk++] = idM;
                surfinfo[sk++] = idP;
                surfinfo[sk++] = (float)( sJk[f]/( 2.0*J[shape->Fmask[f][m]] ) );
                surfinfo[sk++] = (float)sJk[f];
                surfinfo[sk++] = (float)nxk[f];
                surfinfo[sk++] = (float)nyk[f];
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

}
