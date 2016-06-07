#include "Convection2d/Convection2d.h"
#include <mpi.h>

typedef struct face {
    int p1, k1, f1;
    int p2, k2, f2;
    /** max index of vertex on face */
    int va;
    /** min index of vertex on face */
    int vb;
    /** max index of vertex */
    int g;
}face;

/**
 * @brief
 * Comparation of face objects: obj1 and obj2
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 * @return
 * return values：
 * value| description of value
 * ---  |-------
 * -1   | not adjacent
 *  1   | not adjacent
 *  0   | two faces are adjacent
 */
int compare_pairs(const void *obj1, const void *obj2){

    face *e1 = (face*) obj1;
    face *e2 = (face*) obj2;

    int a1 = e1->va, b1 = e1->vb;
    int a2 = e2->va, b2 = e2->vb;

    int va1, vb1, va2, vb2;

    va1 = min(a1, b1);
    vb1 = max(a1, b1);

    va2 = min(a2, b2);
    vb2 = max(a2, b2);

    if(vb1<vb2)
        return -1;
    else if(vb1>vb2)
        return 1;
    else if(va1<va2)
        return -1;
    else if(va1>va2)
        return 1;

    return 0;
}


int pairprocget(const void *obj1){
    face *e1 = (face*) obj1;
    return (e1->p1);
}


int pairnumget(const void *obj1){
    face *e1 = (face*) obj1;
    return (e1->g);
}

void pairnumset(const void *obj1, int g){
    face *e1 = (face*) obj1;
    e1->g = g;
}

void pairmarry(const void *obj1, const void *obj2){
    face *e1 = (face*) obj1;
    face *e2 = (face*) obj2;
    e1->p2 = e2->p1;  e1->k2 = e2->k1;  e1->f2 = e2->f1;
    e2->p2 = e1->p1;  e2->k2 = e1->k1;  e2->f2 = e1->f1;
}

/**
 * @brief 简要说明
 * @details 详细说明
 * @date
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 * @param[in] beginPos 对应区域开始显示的地址
 * @param[in] order order>0: year/month/date;order=0: date/month/year
 * @return
 * return values：
 * name     | type     | description of value
 * -------- |----------|----------------------
 * car_id   | int      |
 * car_info | object   |
 *
 */
void FacePair(Mesh *mesh) {

    int procid = mesh->procid;
    int nprocs = mesh->nprocs;

    int Klocal = mesh->K;
    int Nfaces = mesh->Nfaces;

    int **EToV = mesh->EToV;

#if defined TRI
    const int vnum[3][2] = { {0,1}, {1,2}, {2,0} };
#elif defined QUAD
    const int vnum[4][2] = { {0,1}, {1,2}, {2,3}, {3,0}};
#endif
    int n, k, e, sk;

    face *myfaces = (face*) calloc(Klocal*Nfaces, sizeof(face));

    sk = 0;
    for(k=0;k<Klocal;++k){
        for(e=0;e<Nfaces;++e){
            int n1 = EToV[k][vnum[e][0]];
            int n2 = EToV[k][vnum[e][1]];

            myfaces[sk].p1 = procid; myfaces[sk].k1 = k; myfaces[sk].f1 = e;
            myfaces[sk].p2 = procid; myfaces[sk].k2 = k; myfaces[sk].f2 = e;

            myfaces[sk].va = max(n1,n2);
            myfaces[sk].vb = min(n1,n2);
            myfaces[sk].g  = max(n1,n2); /* marker for sorting into bins */
            ++sk;
        }
    }

    ParallelPairs(myfaces, Klocal*Nfaces, sizeof(face),
                  pairnumget, pairnumset, pairprocget,
                  pairmarry, compare_pairs);

    mesh->Npar = BuildIntVector(nprocs);

    mesh->EToE = BuildIntMatrix(Klocal, Nfaces);
    mesh->EToF = BuildIntMatrix(Klocal, Nfaces);
    mesh->EToP = BuildIntMatrix(Klocal, Nfaces);

    int k1, k2, f1, f2, p1, p2;

    for(n=0;n<Klocal*Nfaces;++n){

        k1 = myfaces[n].k1; f1 = myfaces[n].f1; p1 = myfaces[n].p1;
        k2 = myfaces[n].k2; f2 = myfaces[n].f2; p2 = myfaces[n].p2;

        if(p1!=procid)
            fprintf(stderr, "WARNING WRONG proc\n");

        mesh->EToE[k1][f1] = k2;
        mesh->EToF[k1][f1] = f2;
        mesh->EToP[k1][f1] = p2;

        if(p1!=p2){
            /* increment number of links */
            ++mesh->Npar[p2];
        }
    }

    mesh->parK = (int**) calloc(nprocs, sizeof(int*));
    mesh->parF = (int**) calloc(nprocs, sizeof(int*));
    for(p2=0;p2<nprocs;++p2){
        mesh->parK[p2] = BuildIntVector(mesh->Npar[p2]);
        mesh->parF[p2] = BuildIntVector(mesh->Npar[p2]);
        mesh->Npar[p2] = 0;
        for(n=0;n<Klocal*Nfaces;++n){
            if(myfaces[n].p2==p2 && p2!=procid){
                k1 = myfaces[n].k1, f1 = myfaces[n].f1;
//                k2 = myfaces[n].k2, f2 = myfaces[n].f2;
                mesh->parK[p2][mesh->Npar[p2]  ] = k1;
                mesh->parF[p2][mesh->Npar[p2]++] = f1;
            }
        }
    }

    free(myfaces);
}

#define DSET_NAME_LEN 1024

/**
 * @brief
 * print out mesh connectivity
 *
 * @details
 * each processor generate a file and print the result
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 */
void PrintMeshConnection(Mesh* mesh){

    int n, m, rank, nprocs, ret;
    char filename[DSET_NAME_LEN];

    rank = mesh->procid;
    nprocs = mesh->nprocs;

    ret = snprintf(filename, DSET_NAME_LEN, "%d-%d.txt", rank, nprocs);
    if (ret >= DSET_NAME_LEN) {
        fprintf(stderr, "name too long \n");
        exit(-1);
    }

    FILE *fig = fopen(filename, "w");

    fprintf(fig, "Mesh data: \n");
    fprintf(fig, "\n K = %d\n", mesh->K);
    fprintf(fig, "\n Nv = %d\n", mesh->Nv);
    fprintf(fig, "\n Nverts = %d\n", mesh->Nverts);
    fprintf(fig, "\n Element to element connectivity: \n");
    for (n = 0; n < mesh->K; n ++){
        for( m = 0; m < mesh->Nfaces; m++)
            fprintf(fig, "%d,\t", mesh->EToE[n][m]);
        fprintf(fig, "\n");
    }

    fprintf(fig, "\n Element to face connectivity: \n");
    for (n = 0; n < mesh->K; n ++){
        for( m = 0; m < mesh->Nfaces; m++)
            fprintf(fig, "%d,\t", mesh->EToF[n][m]);
        fprintf(fig, "\n");
    }

    fprintf(fig, "\n Element to processor connectivity = \n");
    for(n=0;n<mesh->K;++n){
        for( m = 0; m < mesh->Nfaces; m++)
            fprintf(fig, "%d,\t", mesh->EToP[n][m]);
        fprintf(fig, "\n");
    }

#if 0
    printf("\n vmapM = \n");
    for(m = 0; m<mesh->K; m++){
        printf("\nElement %d:\n", m);
        for(f1=0;f1<p_Nfaces;++f1) {
            for (n = 0; n < p_Nfp; ++n) {
                int id1 = n + f1 * p_Nfp + m * p_Nfp * p_Nfaces;
                printf("%d, ", mesh->vmapM[id1]);
            }
        }
    }

    printf("\n vmapP = \n");
    for(m = 0; m<mesh->K; m++){
        printf("\nElement %d:\n", m);
        for(f1=0;f1<p_Nfaces;++f1) {
            for (n = 0; n < p_Nfp; ++n) {
                int id1 = n + f1 * p_Nfp + m * p_Nfp * p_Nfaces;
                printf("%d, ", mesh->vmapP[id1]);
            }
        }
    }
#endif

    fclose(fig);
}