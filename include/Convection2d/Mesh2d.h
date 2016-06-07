/* Define physical constant and mesh structure
 *
 * Define physical constant and mesh structure
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 * */

#ifndef MESH2D_H
#define MESH2D_H

/* default degree */
#ifndef p_N
/** default order */
#define p_N 1
#endif

/* default element number in mesh */
#ifndef Ne
#define Ne 80
#endif

#define NODETOL   1e-4

/* element shape constant for triangle */
#ifdef TRI
/** number of nodes in each face */
#define p_Nfp     (p_N+1)
/** number of nodes in elements */
#define p_Np      ((p_N+1)*(p_N+2)/2)
/** number of faces */
#define p_Nfaces  3
#endif

/* element shape constant for quadrilateral */
#ifdef QUAD
/** number of nodes in each face */
#define p_Nfp     (p_N+1)
/** number of nodes in elements */
#define p_Np      ((p_N+1)*(p_N+1))
/** number of faces */
#define p_Nfaces  4
#endif

//#define BSIZE   (16*((p_Np+15)/16))
#define BSIZE p_Np

#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

/** Mesh struct for RKDG methods */
struct RKDG2d {
    /** number of processes */
    int nprocs;
    /** index of this process */
    int procid;

    /** number of vertices per element */
    int Nverts;
    /** number of faces per element */
    int Nfaces;
    /** number of faces per element (3d only) */
    int Nedges;

    /** number of mesh vertices */
    int Nv;
    /** number of mesh elements */
    int K;
    /** element to global vertex list  */
    int **EToV;
    /** element to global face number list */
    int **EToG;
    /** element to global edge number list */
    int **EToS;

    /** element to neighbor element (elements numbered by their proc) */
    int **EToE;
    /** element to neighbor face    (element local number 0,1,2) */
    int **EToF;
    /** element to neighbor proc    (element prc number 0,1,.., nprocs-1) */
    int **EToP;

    /** element to process local vertex list */
    int **localEToV;
    /** number of unique nodes on this process */
    int localNunique;

    /** vector. entry n is 1 if vertex n is on a boundary */
    int *bcflag;

    /** x-coordinates of element vertices */
    double **GX;
    /** y-coordinates of element vertices */
    double **GY;
    /** z-coordinates of element vertices (3d) */
    double **GZ;

    /* high order node info */
    /** node numbers at faces */
    int   **Fmask;
    /* (r,s) coordinates of reference nodes */
    double  *r,   *s,   *t;
    /** local nodal derivative matrices */
    double **Dr, **Ds, **Dt;
    /** local lift matrix */
    double **LIFT;

    /** x-coordinates of element nodes */
    double **x;
    /** y-coordinates of element nodes */
    double **y;
    /** z-coordinates of element nodes (3d) */
    double **z;

    /** volume id of -ve trace of face node */
    int     *vmapM;
    /** volume id of +ve trace of face node */
    int     *vmapP;
    /** # of faces to send recv to each proc */
    int     *Npar;
    /** element of parallel nodes */
    int    **parK;
    /** face of parallel nodes */
    int    **parF;

    /* MPI stuff */

    /** total number of nodes to send */
    int    parNtotalout;
    /** list of nodes to send out */
    int   *parmapOUT;
    /** device list */
    int   *c_parmapOUT;

    /* MPI data buffers */
    float *f_outQ, *f_inQ;

    /* float version of operators */
    float *f_Dr, *f_Ds, *f_Dt;
    float *f_LIFT;

    /* float geometric info */
    /** volume geometric factors */
    float   *vgeo;
    /** surface geometric factors */
    float   *surfinfo;

    /** trouble cell indicator, entry n is 1 if cell is trouble cell */
    int *tcflag;

    /* float field storage (CPU) */
    float  *f_Q, *f_rhsQ, *f_resQ;

    /** flow rate field */
    float  *f_s;

    /* time stepping constants */
    double *rk4a, *rk4b, *rk4c;
};

typedef struct RKDG2d Mesh;

#endif