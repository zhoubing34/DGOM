/**
 * @file
 * Set physical variables
 *
 * @brief
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#include <MultiRegions/mr_mesh.h>
#include "pf_phys.h"

#define DEBUG 0

/* set the surface information of physField object */
static void phys_surfInfo2d(physField *phys);
/* set the volume geometry information `vgeo` */
static void phys_volumeInfo2d(physField *phys);
/* permute the nodal buffers to send/recv */
static void permuteNodalBuffer(physField *phys);
/* permute the cell buffers to send/recv */
static void permuteCellBuffer(physField *phys);

physField* pf_create(int Nfields, parallMesh *mesh){

    physField *phys = (physField *) calloc(1, sizeof(physField));

    phys->dim = mesh->dim;
    phys->mesh = mesh;
    phys->region = mesh->region;
    phys->grid = mesh->grid;
    phys->cell = mesh->cell;

    phys->Nfield = Nfields; ///< number of fields

    const int K = phys->grid->K;
    const int Np = phys->cell->Np;

    /* nodal array */
    phys->f_Q    = (real *) calloc((size_t) K*Np*Nfields, sizeof(real));
    phys->f_rhsQ = (real *) calloc((size_t) K*Np*Nfields, sizeof(real));
    phys->f_resQ = (real *) calloc((size_t) K*Np*Nfields, sizeof(real));
    phys->f_ext  = (real *) calloc((size_t) K*Np*Nfields, sizeof(real)); ///> external data

    /* volume array */
    phys->c_Q = (real *) calloc((size_t) K*Nfields, sizeof(real));

    /* send/recv buffer */
    permuteNodalBuffer(phys);
    permuteCellBuffer(phys);

    /* surface information */
    phys_surfInfo2d(phys);

    /* volume geometry */
    phys_volumeInfo2d(phys);

    /* initial LDG solver */
    phys->viscosity = NULL;

    return phys;
}

void pf_free(physField *phys){

    free(phys->surfinfo);
    free(phys->vgeo);

    free(phys->f_Q);
    free(phys->f_resQ);
    free(phys->f_rhsQ);
    free(phys->f_ext);
    free(phys->c_Q);

    free(phys->f_inQ);
    free(phys->f_outQ);

    free(phys->c_inQ);
    free(phys->c_outQ);

    vector_int_free(phys->nodeIndexOut);
    vector_int_free(phys->cellIndexOut);
}



/**
 * @brief permute the nodal buffers to send/recv
 * @param[in] phys physField object
 */
static void permuteNodalBuffer(physField *phys){

    const int Nfields = phys->Nfield;
    const int parallNodalNum = phys->mesh->parallNodeNum*Nfields;
    int *nodeIndexOut = vector_int_create(parallNodalNum);
    phys->parallNodeNum = parallNodalNum;
    phys->nodeIndexOut = nodeIndexOut;

    size_t sz = (size_t) parallNodalNum;
    phys->f_inQ  = (real *) calloc(sz, sizeof(real));
    phys->f_outQ = (real *) calloc(sz, sizeof(real));

    parallMesh *mesh = phys->mesh;
    int p2,n1,m,fld;
    int nprocs = mesh->nprocs;
    int procid = mesh->procid;
    stdCell *shape = phys->cell;
    int Nfp = shape->Nfp;

    int sp=0, sk=0;
    for(p2=0;p2<nprocs;++p2){
        if(p2!=procid) {
            /* for each received face */
            for (m=0;m<mesh->Npar[p2];++m) {
                for (n1=0;n1<Nfp;++n1) {
                    for (fld=0;fld<Nfields;++fld) {
                        /* sk and node index determine the map relationship */
                        nodeIndexOut[sp++] = Nfields*(mesh->nodeIndexOut[sk])+fld;
                    }
                    sk++;
                }
            }
        }
    }
    return;
}

/**
 * @brief permute the elemental buffers to send/recv
 * @param[in] phys physField object
 */
static void permuteCellBuffer(physField *phys){

    const int Nfields = phys->Nfield;
    const int phys_parallCellNum = phys->mesh->parallCellNum*Nfields;
    int *cellIndexOut = vector_int_create(phys_parallCellNum);
    //int *cellIndexIn = IntVector_create(phys_parallCellNum);

    phys->parallCellNum = phys_parallCellNum;
    //phys->cellIndexIn = cellIndexIn;
    phys->cellIndexOut = cellIndexOut;

    size_t sz = (size_t) phys_parallCellNum;
    phys->c_inQ  = (real *) calloc(sz, sizeof(real));
    phys->c_outQ = (real *) calloc(sz, sizeof(real));

    parallMesh *mesh = phys->mesh;

    int p2, m, fld;
    int nprocs = mesh->nprocs;
    int procid = mesh->procid;
    int sp=0, sk=0;
    for(p2=0;p2<nprocs;++p2){
        if(p2!=procid) {
            /* for each received face */
            for (m=0;m<mesh->Npar[p2];++m) {
                for (fld = 0; fld < Nfields;++fld) {
                    /* sk and node index determine the map relationship */
                    cellIndexOut[sp] = Nfields*(mesh->cellIndexOut[sk]) + fld;
                    //cellIndexIn[sp] = Nfields*(mesh->cellIndexIn[sk]) + fld;
                    sp++;
                }
                sk++;
            }
        }
    }
}

/**
 * @brief set the volume geometry information `vgeo`
 * @param[in] phys physField object
 */
static void phys_volumeInfo2d(physField *phys){

    multiReg *region = phys->region;
    const int K = phys->grid->K;
    const int Np = phys->cell->Np;

    int Nvgeo = 4;
    size_t sz = (size_t) K*Np*Nvgeo;
    real *vgeo = (real*) calloc(sz, sizeof(real));

    phys->Nvgeo = Nvgeo;
    phys->vgeo = vgeo;

    int k,n,sk=0;
    for(k=0;k<K;k++){
        double *drdx = region->drdx[k];
        double *drdy = region->drdy[k];
        double *dsdx = region->dsdx[k];
        double *dsdy = region->dsdy[k];
        for(n=0;n<Np;n++){
            vgeo[sk++] = drdx[n];
            vgeo[sk++] = drdy[n];
            vgeo[sk++] = dsdx[n];
            vgeo[sk++] = dsdy[n];
        }
    }
}

/**
 * @brief set the surface information of physField object
 * @param[in,out] phys physField object
 */
static void phys_surfInfo2d(physField *phys){

    int k,f,m,sk = 0;
    parallMesh *mesh = phys->mesh;
    stdCell *shape = mesh->cell;
    multiReg *region = phys->region;

    const int Nfields = phys->Nfield;
    const int K = phys->grid->K;
    const int Np = shape->Np;
    const int Nfaces = shape->Nfaces;
    const int Nfp = shape->Nfp;
    int **Fmask = shape->Fmask;

    // allocation of surfinfo
    int Nsurfinfo = 6;
    size_t sz = (size_t) K*Nfp*Nfaces*Nsurfinfo;
    real *surfinfo = (real*) calloc(sz, sizeof(real));
    phys->Nsurfinfo = Nsurfinfo;
    phys->surfinfo = surfinfo;

    if( mesh->EToBS == NULL ){
        fprintf(stderr, "%s(%s): %d\n "
                "mesh->EToBS is not allocated, please set surface types\n",
                __FUNCTION__,__FILE__, __LINE__);
    }
    for(k=0;k<K;++k){

        double *sJk = region->sJ[k];
        double *nxk = region->nx[k];
        double *nyk = region->ny[k];

        double *J = region->J[k];

        for(f=0;f<Nfaces;++f){
            int bcType = mesh->EToBS[k][f];
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
                    idP = Nfields*(-1-idP);
                }

                surfinfo[sk++] = idM; ///< local node index
                surfinfo[sk++] = idP; ///< adjacent node index
                surfinfo[sk++] = (real)(sJk[f]/(J[ Fmask[f][m] ])); ///< face scale
                surfinfo[sk++] = bcType; ///< boundary type indicator
                surfinfo[sk++] = (real)nxk[f]; ///< outward vector, nx
                surfinfo[sk++] = (real)nyk[f]; ///< outward vector, ny
            }
        }
    }
}
