#include "dg_mesh.h"
#include "dg_mesh_connect.h"
#include "dg_mesh_fetch_buffer.h"

#define DEBUG 0

/* sort numbers from small to large */
static int dg_mesh_face_cmp(const void *obj1, const void *obj2);

typedef struct face2d{
    long int tmp;
    long int k1, k2; ///< local and adjacent cell index
    int f1, f2; ///< local and adjacent face index
}face2d;


void dg_mesh_init_node_fetch_buffer(dg_mesh *mesh){
    dg_cell *cell = mesh->cell;
    const int Np = dg_cell_Np(cell);
    int K = dg_grid_K(mesh->grid);
    const int Nparf = dg_mesh_Nparf(mesh);
    const int nprocs = dg_mesh_nprocs(mesh);
    const int procid = dg_mesh_procid(mesh);

    int *parface = mesh->parface;
    // total number of parallel face points
    int f, Nparn = 0;
    for(f=0;f<Nparf;f++){
        Nparn += dg_cell_Nfp(cell)[parface[f]];
    }
    // adjacent node id
    int *n1 = vector_int_create(Nparn);
    int *n2 = vector_int_create(Nparn);
    double *xM = vector_double_create(Nparn);
    double *yM = vector_double_create(Nparn);
    double *xP = vector_double_create(Nparn);
    double *yP = vector_double_create(Nparn);
    int *n_recv = vector_int_create(Nparn);

    int **Fmask = dg_cell_Fmask(cell);
    const int *parfN = mesh->parfaceNum;
    double **x = mesh->region->x;
    double **y = mesh->region->y;
    int Kprocs[nprocs], Ntotal;
    // get number of cells in each process
    MPI_Allgather(&K, 1, MPI_INT, Kprocs, 1, MPI_INT, MPI_COMM_WORLD);

    /* assignment of vampM */
    int k,n,m,t,p,sk=0;
    int *parcell = mesh->parcell;
    for(n=0;n<Nparf;n++){
        k = parcell[n];
        f = parface[n];
        int Nfp = dg_cell_Nfp(cell)[f];
        for(m=0;m<Nfp;m++){
            int ind = k*Np + Fmask[f][m];
            n1[sk] = ind;
            xM[sk] = x[0][ind];
            yM[sk] = y[0][ind];
            sk++;
        }
    }

    /* send to adjacent process */
    MPI_Request vmap_send_requests[nprocs], vmap_recv_requests[nprocs];
    MPI_Request x_send_requests[nprocs], x_recv_requests[nprocs];
    MPI_Request y_send_requests[nprocs], y_recv_requests[nprocs];
    int Nmess=0, st=0, parnodeNum[nprocs]; // number of points adjacent to each process
    sk=0;

    for(p=0;p<nprocs;++p){
        parnodeNum[p] = 0; // initialize
        if(p!=procid){
            int Nout = 0; // # of points send to process p
            for(f=0;f<(parfN[p]);f++){
                Nout += dg_cell_Nfp(cell)[parface[st++]];
            }
            //const int Nout = mesh->parfaceNum[p]*Nfp;
            parnodeNum[p] = Nout;
            if(Nout){
                /* symmetric communications (different ordering) */
                MPI_Isend(n1+sk, Nout, MPI_INT, p, 2333+p, MPI_COMM_WORLD, vmap_send_requests +Nmess);
                MPI_Isend(xM+sk, Nout, MPI_DOUBLE, p, 4333+p, MPI_COMM_WORLD, x_send_requests +Nmess);
                MPI_Isend(yM+sk, Nout, MPI_DOUBLE, p, 5333+p, MPI_COMM_WORLD, y_send_requests +Nmess);
                MPI_Irecv(n_recv+sk,  Nout, MPI_INT, p, 2333+procid, MPI_COMM_WORLD,  vmap_recv_requests +Nmess);
                MPI_Irecv(xP+sk,  Nout, MPI_DOUBLE, p, 4333+procid, MPI_COMM_WORLD,  x_recv_requests +Nmess);
                MPI_Irecv(yP+sk,  Nout, MPI_DOUBLE, p, 5333+procid, MPI_COMM_WORLD,  y_recv_requests +Nmess);
                sk+=Nout;
                ++Nmess;
            }
        }
    }
    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, vmap_recv_requests, instatus);
    MPI_Waitall(Nmess, vmap_send_requests, instatus);
    MPI_Waitall(Nmess, x_recv_requests, instatus);
    MPI_Waitall(Nmess, x_send_requests, instatus);
    MPI_Waitall(Nmess, y_recv_requests, instatus);
    MPI_Waitall(Nmess, y_send_requests, instatus);

    /* set n2, n2 the node index from adjacent process */
    for(n=0;n<Nparf;n++){
        int Nfp = dg_cell_Nfp(cell)[parface[n]];
        for(m=0;m<Nfp;m++){
            sk = n*Nfp+m;
            double x1 = xM[sk];
            double y1 = yM[sk];
            for(t=0;t<Nfp;t++){
                st = n*Nfp+t;
                double x2 = xP[st];
                double y2 = yP[st];
                double d12 = ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
                if(d12 <= EPS) n2[sk] = n_recv[st];
            }
        }
    }

    /* sort the node sequence */
    face2d *my_node = (face2d*) calloc(Nparn, sizeof(face2d));
    sk = 0;
    for(p=0;p<nprocs;p++){
        Ntotal = ((procid) < (p) ? (K*Np) : (Kprocs[p]*Np));
        for(n=0;n<parnodeNum[p];n++){
            long int k1 = (long int)n1[sk];
            long int k2 = (long int)n2[sk];
            my_node[sk].k1 = k1;
            my_node[sk].k2 = k2;
            my_node[sk].tmp = ((procid) < (p) ? (k2*Ntotal+k1) : (k1*Ntotal+k2));
            sk++;
        }
    }

    int stride = 0;
    for(p=0;p<nprocs;p++){
        if(parfN[p]){ // sort the parallel faces with tmp
            qsort(my_node+stride, (size_t)parnodeNum[p], sizeof(face2d), dg_mesh_face_cmp);
        }
        stride += parnodeNum[p];
    }
    /* assignment */
    mesh->TotalParNode = Nparn;
    mesh->parnodeNum = (int *) calloc(nprocs, sizeof(int));
    mesh->parnode = (int *) calloc(Nparn, sizeof(int));

    for(p=0;p<nprocs;p++){
        mesh->parnodeNum[p] = parnodeNum[p];
    }
    for(n=0;n<Nparn;n++){
        mesh->parnode[n] = (int)my_node[n].k1;
    }
    /* free */
    vector_double_free(xM);
    vector_double_free(xP);
    vector_double_free(yM);
    vector_double_free(yP);
    vector_int_free(n1);
    vector_int_free(n2);
    vector_int_free(n_recv);
    free(my_node);
    return;
}

/**
 * @brief
 * @param mesh
 */
void dg_mesh_init_cell_fetch_buffer(dg_mesh *mesh){

    const dg_grid *grid = mesh->grid;
    const int nprocs = mesh->nprocs;
    const int procid = mesh->procid;
    const int Nfaces = dg_cell_Nfaces(mesh->cell);

    int Klocal = dg_grid_K(mesh->grid);
    /* number of faces adjacent to other process */
    int Nparf=0, Parf[nprocs];

    int k,f,n,p2;
    for(p2=0;p2<nprocs;p2++){
        Parf[p2] = 0;
        for(k=0;k<Klocal;k++){
            for(f=0;f<Nfaces;f++){
                if( (grid->EToP[k][f] != procid) & (grid->EToP[k][f] == p2) ){
                    /* increment number of links */
                    Parf[p2] += 1;
                    Nparf++;
                }
            }
        }
    }

    int Kprocs[nprocs], Ktotal;
    // get number of cells in each process
    MPI_Allgather(&Klocal, 1, MPI_INT, Kprocs, 1, MPI_INT, MPI_COMM_WORLD);
    face2d *my_face = (face2d *)calloc(Nparf, sizeof(face2d));

    int sk=0;
    for(p2=0;p2<nprocs;p2++){
        if(Parf[p2]){
            // cell number for process with small index
            Ktotal = ((procid) < (p2) ? (Klocal) : (Kprocs[p2]));
            for(k=0;k<Klocal;k++){
                for(f=0;f<Nfaces;f++){
                    if(grid->EToP[k][f]==p2){
                        long int k1 = k;//
                        long int k2 = grid->EToE[k][f];// + Kstart2;
                        my_face[sk].k1 = k1;
                        my_face[sk].k2 = k2;
                        my_face[sk].f1 = f;
                        my_face[sk].f2 = grid->EToF[k][f];
                        // tmp is the indicator for each face pairs
                        // k2 * Ktotal + k1
                        // k1 is the cell index in process 1 (smaller process id)
                        // k2 and Ktotal is the cell index and total cell number in process 2 (greater process id)
                        my_face[sk].tmp = ((procid) < (p2) ? (k2*Ktotal+k1) : (k1*Ktotal+k2));
                        sk++;
                    }
                }
            }
        }
    }

    int stride = 0;
    for(p2=0;p2<nprocs;p2++){
        if(Parf[p2]){ // sort the parallel faces with tmp
            qsort(my_face+stride, (size_t)Parf[p2], sizeof(face2d), dg_mesh_face_cmp);
        }
        stride += Parf[p2];
    }

    /* assignment */
    mesh->TotalParFace = Nparf;
    mesh->parfaceNum = vector_int_create(nprocs);
    mesh->parcell = vector_int_create(Nparf);
    mesh->parface = vector_int_create(Nparf);

    for(n=0;n<Nparf;n++){
        mesh->parcell[n] = (int)my_face[n].k1;
        mesh->parface[n] = my_face[n].f1;
    }
    for(p2=0;p2<nprocs;p2++){
        mesh->parfaceNum[p2] = Parf[p2];
    }

    free(my_face);
    return;
}

/* sort numbers from small to large */
static int dg_mesh_face_cmp(const void *obj1, const void *obj2){
    face2d *e1 = (face2d*) obj1;
    face2d *e2 = (face2d*) obj2;
#if DEBUG
    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    printf("procid=%d, %s (%d), e1: k1=%ld, k2=%ld, f1=%d, f2=%d, tmp=%ld\n",
           procid, __FUNCTION__, __LINE__, e1->k1, e1->k2, e1->f1, e1->f2, e1->tmp);
    printf("procid=%d, %s (%d), e2: k1=%ld, k2=%ld, f1=%d, f2=%d, tmp=%ld\n",
           procid, __FUNCTION__, __LINE__, e2->k1, e2->k2, e2->f1, e2->f2, e2->tmp);
#endif
    if ( (e1->tmp) > (e2->tmp) ) {return 1;}
    else if( (e1->tmp) < (e2->tmp) ) {return -1;}
    else if( (e1->tmp) == (e2->tmp) ) {
        int procid;
        MPI_Comm_rank(MPI_COMM_WORLD, &procid);
        fprintf(stderr, "procid=%d, %s (%d)\nCell pairs %ld-%ld and %ld-%ld, the indicator should not be equal.\n",
                procid, __FUNCTION__, __LINE__, e1->k1, e1->k2, e2->k1, e2->k2);
        exit(-1);
    }
    return 0;
}