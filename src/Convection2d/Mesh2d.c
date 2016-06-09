#include "Convection2d/Convection2d.h"
#include <mpi.h>

/**
 * @brief
 * Generate local uniform triangle mesh for each processor
 *
 * @details
 * The whole domain is [-1, 1]x[-1, 1], with each edge partition to *Ne* 
 * elements. The vertex list is arranged from upper left to lower right
 * directions. The whole mesh is generate with newEToV for element to 
 * vertex list. Each process reads part of mesh into mesh->EToV and its 
 * corresponding vertex coordinate mesh->GX and mesh->GY. The uniform 
 * triangle element is separated from the square, and the separation is 
 * controlled by UP_RIGHT macro.
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 * @return
 * return values：
 * name     | type      | description of value
 * -------- |---------- |----------------------
 * mesh     | Mesh*     | mesh object
 *
 */
Mesh* ReadTriMesh(){

    /* mesh object and allocation */
    Mesh *mesh = (Mesh*) calloc(1, sizeof(Mesh));
    /* index of vertex at upper and lower layer */
    int upperNode[2], lowerNode[2];
    int irow, icol, ie, p, K;
    /** number of elements in local processor */
    int Klocal;
    /** start index of element in all element */
    int Kstart = 0;
    /** number of elements in each processor */
    int * Kprocs;

    /* decide on parition */
    int procid, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if(!procid) printf("Root: Entering ReadTriMesh\n");

    mesh->procid = procid;
    mesh->nprocs = nprocs;

    mesh->Nverts = 3; /* assume triangles */
    mesh->Nfaces = 3; /* assume triangles */

    Kprocs = (int *) calloc(nprocs, sizeof(int));
    K = Ne* Ne *2; /* total # of elements */
    Klocal = (int)( (double)K/ (double)nprocs );

    /* number of elements in each process */
    for(p = 0; p < nprocs - 1; ++p){
        Kprocs[p] = Klocal;
    }
    /* the rest elements is distributed to the last process: p = nprocs - 1  */
    Kprocs[p] = K + Klocal - nprocs*Klocal;

    /* number of element in local process */
    Klocal = Kprocs[procid];

    /* start index of element */
    for(p=0;p<procid;++p)
        Kstart += Kprocs[p];

    /* # of all vertex */
    mesh->Nv = (Ne + 1)*(Ne + 1);

    /* allocate coordinate of all vertex */
    double *VX = BuildVector(mesh->Nv);
    double *VY = BuildVector(mesh->Nv);

    /* element to vertex list for all elements */
    int **newEToV = BuildIntMatrix(K, mesh->Nverts);
    /* EToV of local mesh */
    mesh->EToV = BuildIntMatrix(Klocal, mesh->Nverts);

    for (irow = 0; irow < Ne; irow++){
        for (ie = 0; ie < Ne; ie++){
            /* Assignment of vertex index, which start from 0 */
            upperNode[0] = irow*(Ne + 1)+ie;
            upperNode[1] = irow*(Ne + 1)+ ie + 1;
            lowerNode[0] = upperNode[0] + Ne + 1;
            lowerNode[1] = upperNode[1] + Ne + 1;

            /* Assignment of EToV */
#define UP_RIGHT
#ifdef UP_RIGHT
            newEToV[irow * Ne * 2 + ie][0] = lowerNode[0];
            newEToV[irow * Ne * 2 + ie][1] = upperNode[1];
            newEToV[irow * Ne * 2 + ie][2] = upperNode[0];

            newEToV[irow * Ne * 2 + Ne + ie][0] = lowerNode[0];
            newEToV[irow * Ne * 2 + Ne + ie][1] = lowerNode[1];
            newEToV[irow * Ne * 2 + Ne + ie][2] = upperNode[1];
#else // up-left
            newEToV[irow * Ne * 2 + ie][0] = lowerNode[1];
            newEToV[irow * Ne * 2 + ie][1] = upperNode[1];
            newEToV[irow * Ne * 2 + ie][2] = upperNode[0];

            newEToV[irow * Ne * 2 + Ne + ie][0] = lowerNode[0];
            newEToV[irow * Ne * 2 + Ne + ie][1] = lowerNode[1];
            newEToV[irow * Ne * 2 + Ne + ie][2] = upperNode[0];
#endif
        }
    }

    /* Assignment of coordinate for vertex */
    int index; /* vertex index */
    for (irow = 0; irow<Ne+1; irow++){
        for (icol = 0; icol<Ne+1; icol++){
            index = irow*(Ne+1) + icol;
            VX[index] = (double)icol/Ne *2.0 -1;
            VY[index] = (double)irow/Ne * (-2.0) + 1;
        }
    }

    mesh->GX = BuildMatrix(Klocal, mesh->Nverts);
    mesh->GY = BuildMatrix(Klocal, mesh->Nverts);

    /* Assignment of local mesh information */
    int sk = 0;
    for (ie = 0; ie < K; ie++){
        if(ie>=Kstart && ie<Kstart+Klocal) {
            mesh->EToV[sk][0] = newEToV[ie][0];
            mesh->EToV[sk][1] = newEToV[ie][1];
            mesh->EToV[sk][2] = newEToV[ie][2];

            mesh->GX[sk][0] = VX[mesh->EToV[sk][0]];
            mesh->GX[sk][1] = VX[mesh->EToV[sk][1]];
            mesh->GX[sk][2] = VX[mesh->EToV[sk][2]];

            mesh->GY[sk][0] = VY[mesh->EToV[sk][0]];
            mesh->GY[sk][1] = VY[mesh->EToV[sk][1]];
            mesh->GY[sk][2] = VY[mesh->EToV[sk][2]];
            ++sk;

        }
    }

    mesh->K = Klocal;

    DestroyVector(VX);
    DestroyVector(VY);
    DestroyIntMatrix(newEToV);

    free(Kprocs);

    if(!procid) printf("Root: Leaving ReadTriMesh\n");
    return mesh;
}

/**
 * @brief
 * Generate uniform quadrilateral mesh
 *
 * @details
 * The whole domain is [-1, 1]x[-1, 1], with each edge partition to *Ne* 
 * elements. The vertex list is arranged from upper left to lower right
 * directions. The whole mesh is generate with newEToV for element to 
 * vertex list. Each processor reads part of mesh into mesh->EToV and 
 * its corresponding vertex coordinate mesh->GX and mesh->GY.
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 * @return
 * return values：
 * name     | type      | description of value
 * -------- |---------- |----------------------
 * mesh     | Mesh*     | mesh object
 *
 */
Mesh* ReadQuadMesh(){
    /* mesh object and allocation */
    Mesh *mesh = (Mesh*) calloc(1, sizeof(Mesh));
    /* index of vertex at upper and lower layer */
    int upperNode[2], lowerNode[2];
    int irow, icol, ie, p, K;
    /** number of elements in local processor */
    int Klocal;
    /** start index of element in all element */
    int Kstart = 0;
    /** number of elements in each processor */
    int * Kprocs;

    /* decide on parition */
    int procid, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if(!procid) printf("Root: Entering ReadTriMesh\n");

    mesh->procid = procid;
    mesh->nprocs = nprocs;

    mesh->Nverts = 4; /* assume triangles */
    mesh->Nfaces = 4; /* assume triangles */

    Kprocs = (int *) calloc(nprocs, sizeof(int));
    K = Ne* Ne; /* total # of elements */
    Klocal = (int)( (double)K/ (double)nprocs );

    /* number of elements in each process */
    for(p = 0; p < nprocs - 1; ++p){
        Kprocs[p] = Klocal;
    }
    /* the rest elements is distributed to the last process: p = nprocs - 1  */
    Kprocs[p] = K + Klocal - nprocs*Klocal;

    /* number of element in local process */
    Klocal = Kprocs[procid];

    /* start index of element */
    for(p=0;p<procid;++p)
        Kstart += Kprocs[p];

    /* # of all vertex */
    mesh->Nv = (Ne + 1)*(Ne + 1);

    /* allocate coordinate of all vertex */
    double *VX = BuildVector(mesh->Nv);
    double *VY = BuildVector(mesh->Nv);

    /* element to vertex list for all elements */
    int **newEToV = BuildIntMatrix(K, mesh->Nverts);
    /* EToV of local mesh */
    mesh->EToV = BuildIntMatrix(Klocal, mesh->Nverts);

    for (irow = 0; irow < Ne; irow++){
        for (ie = 0; ie < Ne; ie++){
            /* Assignment of vertex index, which start from 0 */
            upperNode[0] = irow*(Ne + 1)+ie;
            upperNode[1] = irow*(Ne + 1)+ ie + 1;
            lowerNode[0] = upperNode[0] + Ne + 1;
            lowerNode[1] = upperNode[1] + Ne + 1;

            /* Assignment of EToV,
             * each element start from the left lower vertex */
            newEToV[irow * Ne + ie][0] = lowerNode[0];
            newEToV[irow * Ne + ie][1] = lowerNode[1];
            newEToV[irow * Ne + ie][2] = upperNode[1];
            newEToV[irow * Ne + ie][3] = upperNode[0];
        }
    }

    /* Assignment of coordinate for vertex */
    int index; /* vertex index */
    for (irow = 0; irow<Ne+1; irow++){
        for (icol = 0; icol<Ne+1; icol++){
            index = irow*(Ne+1) + icol;
            VX[index] = (double)icol/Ne *2.0 -1;
            VY[index] = (double)irow/Ne * (-2.0) + 1;
        }
    }

    mesh->GX = BuildMatrix(Klocal, mesh->Nverts);
    mesh->GY = BuildMatrix(Klocal, mesh->Nverts);

    /* Assignment of local mesh information */
    int sk = 0;
    for (ie = 0; ie < K; ie++){
        if(ie>=Kstart && ie<Kstart+Klocal) {
            mesh->EToV[sk][0] = newEToV[ie][0];
            mesh->EToV[sk][1] = newEToV[ie][1];
            mesh->EToV[sk][2] = newEToV[ie][2];
            mesh->EToV[sk][3] = newEToV[ie][3];

            mesh->GX[sk][0] = VX[mesh->EToV[sk][0]];
            mesh->GX[sk][1] = VX[mesh->EToV[sk][1]];
            mesh->GX[sk][2] = VX[mesh->EToV[sk][2]];
            mesh->GX[sk][3] = VX[mesh->EToV[sk][3]];

            mesh->GY[sk][0] = VY[mesh->EToV[sk][0]];
            mesh->GY[sk][1] = VY[mesh->EToV[sk][1]];
            mesh->GY[sk][2] = VY[mesh->EToV[sk][2]];
            mesh->GY[sk][3] = VY[mesh->EToV[sk][3]];
            ++sk;
        }
    }

    mesh->K = Klocal;

    /* deallocate memory */
    DestroyVector(VX);
    DestroyVector(VY);
    DestroyIntMatrix(newEToV);

    free(Kprocs);
    if(!procid) printf("Root: Leaving ReadTriMesh\n");
    return mesh;
}


/**
 * @brief
 * Get the geometric factor of kth element
 *
 * @details
 * Get the geometric factor in element *k*, the derived geometric factor 
 * on each node is calculated as
 * \f[ \left. \frac{\partial x}{\partial r} \right|_{\mathbf{r}_i} = 
 * \left. \frac{\partial \mathbf{\varphi}^T_j}{\partial r} \right|_{\mathbf{r}_i} \mathbf{x}_j \f]
 * where \f$ \mathbf{\varphi} \f$ is the vector of Lagrangian basis functions. The same as 
 * \f$ \frac{\partial x}{\partial s}, \frac{\partial y}{\partial r} \frac{\partial y}{\partial s} \f$
 *
 * The Jacobian matrix is to transfer the differentials of \f$ \{x,y\}\f$ to 
 * \f$ \{r,s \}\f$. Using the chain rule:
 * \f[ \begin{bmatrix} dx \cr dy \end{bmatrix} =
 * \begin{bmatrix} 
 * \frac{\partial x}{\partial r} & \frac{\partial x}{\partial s} \cr
 * \frac{\partial y}{\partial r} & \frac{\partial y}{\partial s} \cr
 * \end{bmatrix}
 * \begin{bmatrix} dr \cr ds \end{bmatrix} = 
 * \mathbf{J}^T \begin{bmatrix} dr \cr ds \end{bmatrix} \f]
 *
 * The corresponding Jacobian and inverse Jacobian is 
 * \f[ \mathbf{J} = \frac{\partial (x, y)}{\partial (r, s)} =
 * \begin{bmatrix} 
 * \frac{\partial x}{\partial r} & \frac{\partial x}{\partial s} \cr
 * \frac{\partial y}{\partial r} & \frac{\partial y}{\partial s} \cr
 * \end{bmatrix} = \begin{bmatrix} J_{11} & J_{12} \cr
 * J_{21} & J_{22} \end{bmatrix}, \quad 
 * \mathbf{J}^{-1} = \frac{\partial (r, s)}{\partial (x,y)} =
 * \begin{bmatrix} 
 * \frac{\partial r}{\partial x} & \frac{\partial s}{\partial x} \cr
 * \frac{\partial r}{\partial y} & \frac{\partial s}{\partial y} \cr
 * \end{bmatrix} = \frac{1}{J} \begin{bmatrix} J_{22} & -J_{12} \cr
 * -J_{21} & J_{11} \end{bmatrix}, \quad \f]
 * where \f$ J = det(\mathbf{J}) = J_{11}J_{22} - J_{12}J_{21}\f$
 * 
 * @author 
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 * @param[in] mesh The Mesh object
 * @param[in] k The element index
 *
 * @return  
 * return values：
 * name     | type      | description of value
 * -------- |---------- |----------------------
 * drdx     | double*   | the drdx on each node
 * drdy     | double*   | the drdy on each node
 * dsdx     | double*   | the dsdx on each node
 * dsdy     | double*   | the dsdy on each node
 * J        | double*   | the Jacobi on each node
 *
 */
void GeometricFactors(Mesh *mesh, int k,
                      double *drdx, double *dsdx, double
                      *drdy, double *dsdy, double *J){

    double *dxdr, *dxds;
    double *dydr, *dyds;
    const double *x = mesh->x[k];
    const double *y = mesh->y[k];
    int n, m;

    dxdr = (double *)calloc(p_Np, sizeof(double));
    dxds = (double *)calloc(p_Np, sizeof(double));
    dydr = (double *)calloc(p_Np, sizeof(double));
    dyds = (double *)calloc(p_Np, sizeof(double));

    for (n = 0; n<p_Np; n++){
        for (m = 0; m<p_Np; m++){
            dxdr[n] += mesh->Dr[n][m]*x[m];
            dxds[n] += mesh->Ds[n][m]*x[m];
            dydr[n] += mesh->Dr[n][m]*y[m];
            dyds[n] += mesh->Ds[n][m]*y[m];
        }
        /* Jacobian of coordinate mapping */
        J[n] = -dxds[n]*dydr[n] + dxdr[n]*dyds[n];

        if(J[n]<0)
            printf("warning: J = %lg\n", J[n]);

        /* inverted Jacobian matrix for coordinate mapping */
        drdx[n] =  dyds[n]/(J[n]);
        dsdx[n] = -dydr[n]/(J[n]);
        drdy[n] = -dxds[n]/(J[n]);
        dsdy[n] =  dxdr[n]/(J[n]);
    }

    free(dxdr);
    free(dxds);
    free(dydr);
    free(dyds);
}

/**
 * @brief
 * Get normal vector and Jacobian coefficient of each faces in kth element
 *
 * @details
 * The outward normal vector \f$ \vec{n} = \left(n_x, n_y \right) \f$ is perpendicular
 * to each side \f$ \vec{r} = \left(\Delta x, \Delta y \right) \f$, thus the normal vector
 * can be obtained as
 * \f[ \vec{n} = \frac{1}{s} \left(\Delta y, -\Delta x \right) \f]
 * where s is the length of each side. The Jacobian transfer coefficient is 
 * \f$ J_s = s/2 \f$, where the length of each sides in standard quadrilateral 
 * element is 2.
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 */
void Normals(Mesh *mesh, int k, double *nx, double *ny, double *sJ){
    int f;
    double x1, x2, y1, y2;

    /*  */
    for(f=0; f<(mesh->Nfaces -1); f++){
        x1 = mesh->GX[k][f]; x2 = mesh->GX[k][f+1];
        y1 = mesh->GY[k][f]; y2 = mesh->GY[k][f+1];
        nx[f] =  (y2-y1);
        ny[f] = -(x2-x1);
    }
    /* the last face */
    x1 = mesh->GX[k][f]; x2 = mesh->GX[k][0];
    y1 = mesh->GY[k][f]; y2 = mesh->GY[k][0];
    nx[f] =  (y2-y1);
    ny[f] = -(x2-x1);

    /* normalize the outward vector */
    for(f=0;f<mesh->Nfaces;++f){
        sJ[f] = sqrt(nx[f]*nx[f]+ny[f]*ny[f]);
        nx[f] /= sJ[f];
        ny[f] /= sJ[f];
        sJ[f] /= 2.;
    }
}


/**
 * @brief
 * Print the mesh information of each process in files
 *
 * @details
 *
 * @author 
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 */
#define DSET_NAME_LEN 1024
void PrintMesh ( Mesh *mesh ){
    int n, m, rank, nprocs, ret, sk, sf;
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
    fprintf(fig, "\n Vertex coordinates X: \n");
    for (n = 0; n < mesh->K; n ++){
        for( m = 0; m < mesh->Nverts; m++)
            fprintf(fig, "%f,\t", mesh->GX[n][m]);
        fprintf(fig, "\n");
    }

    fprintf(fig, "\n Vertex coordinates Y: \n");
    for (n = 0; n < mesh->K; n ++){
        for( m = 0; m < mesh->Nverts; m++)
            fprintf(fig, "%f,\t", mesh->GY[n][m]);
        fprintf(fig, "\n");
    }

#if 0
    fprintf(fig, "\n Node coordinate = \n");
    for (n = 0; n < mesh->K; n ++){
        for( m = 0; m < p_Np; m++)
            fprintf(fig, "[%f, %f],\t", mesh->x[n][m], mesh->y[n][m]);
        fprintf(fig, "\n");
    }

#endif

#if 0

    sk = 0;
    fprintf(fig, "\n Node vgeo = \n");
    for (n = 0; n < mesh->K; n ++){
        fprintf(fig, "Ele %d:\n", n);
        for( m = 0; m < p_Np; m++)
            fprintf(fig, "drdx = %f, drdy = %f, dsdx = %f, dsdy = %f\n",
                    mesh->vgeo[sk++], mesh->vgeo[sk++], mesh->vgeo[sk++], mesh->vgeo[sk++]);
        fprintf(fig, "\n");
    }

    sf = 2;
    for (n = 0; n < mesh->K; n ++){
        fprintf(fig, "Ele %d:\n", n);
        for( m = 0; m < p_Nfp*p_Nfaces; m++) {
            fprintf(fig, "Js/J = %f\t", mesh->surfinfo[sf]);
            sf+=6;
        }
        fprintf(fig, "\n");
    }
#endif

    fclose(fig);

}

