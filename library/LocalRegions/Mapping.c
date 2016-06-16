/**
 * \file
 * Mapping the physical elements to standard elements
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "LocalRegions.h"


/**
 * @brief
 * Mapping nodes from standard triangle element to physical element
 * @details
 * \f[ \mathbf{x} = -\frac{r+s}{2}\mathbf{x}_1 + \frac{r+1}{2}\mathbf{x}_2 +
 * \frac{s+1}{2}\mathbf{x}_3 \f]
 *
 * while the coefficient \f$ l_1(r,s) = -\frac{r+s}{2},
 * \quad l_2(r,s) = \frac{r+1}{2},
 * \quad l_3(r,s) = \frac{s+1}{2}\f$ are the linear Lagrangian functions, which
 * fulfill the Lagrange intepolation property
 * \f$ l_i( \mathbf{r}_j ) = \delta_{ij} \f$ in standard element.
 *
 * @param [int]     Np number of points
 * @param [double]  r[Np] coordinate in standard element
 * @param [double]  s[Np] coordinate in standard element
 * @param [double]  GX[3] vertex coordinate
 * @param [double]  GY[3] vertex coordinate
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * x | double[Np]   |  coordinate in physical element
 * y | double[Np]   |  coordinate in physical element
 *
 */

void MapTriCoor(int Np, double *r, double *s,
                double *GX, double *GY,
                double *x, double *y){
    int n;
    double ri, si;
    for (n = 0; n < Np; ++n) {
        ri = r[n]; si = s[n];
        x[n] = 0.5 * (-GX[0] * (ri + si) + GX[1] * (1. + ri) + GX[2] * (1. + si));
        y[n] = 0.5 * (-GY[0] * (ri + si) + GY[1] * (1. + ri) + GY[2] * (1. + si));
    }
}

/**
 * @brief
 * Mapping nodes from standard quadrilateral element to physical element
 * @details
 * \f[ \mathbf{x} = \frac{1-r}{2}\frac{1-s}{2}\mathbf{x}_1 +
 * \frac{1+r}{2}\frac{1-s}{2}\mathbf{x}_2 +
 * \frac{1-r}{2}\frac{1+s}{2}\mathbf{x}_3 +
 * \frac{1+r}{2}\frac{1+s}{2}\mathbf{x}_4 \f]
 *
 * while the coefficient \f$ l_1(r,s) = \frac{1-r}{2}\frac{1-s}{2},
 * \quad l_2(r,s) = \frac{1+r}{2}\frac{1-s}{2},
 * \quad l_3(r,s) = \frac{1-r}{2}\frac{1+s}{2},
 * \quad l_4(r,s) = \frac{1+r}{2}\frac{1+s}{2}\f$ are the
 * linear Lagrangian functions, which fulfill the Lagrange intepolation property
 * \f$ l_i( \mathbf{r}_j ) = \delta_{ij} \f$ in standard element.
 *
 * @param [int]     Np number of points
 * @param [double]  r[Np] coordinate in standard element
 * @param [double]  s[Np] coordinate in standard element
 * @param [double]  GX[3] vertex coordinate
 * @param [double]  GY[3] vertex coordinate
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * x | double[Np]   |  coordinate in physical element
 * y | double[Np]   |  coordinate in physical element
 *
 */
void MapQuadCoor(int Np, double *r, double *s,
                 double *GX, double *GY,
                 double *x, double *y){
    int n;
    double ri, si;
    for (n = 0; n < Np; ++n) {
        ri = r[n];
        si = s[n];
        x[n] = 0.25 * (GX[0] * (1. - ri) * (1. - si) + GX[1] * (1. + ri) * (1. - si)
                       + GX[2] * (1. + ri) * (1. + si) + GX[3] * (1. - ri)*(1. + si));
        y[n] = 0.25 * (GY[0] * (1. - ri) * (1. - si) + GY[1] * (1. + ri) * (1. - si)
                       + GY[2] * (1. + ri) * (1. + si) + GY[3] * (1. - ri)*(1. + si));
    }
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
 * @param [int]      Nv number of vertex
 * @param [double*]  GX[Nv] vertex coordinate
 * @param [double*]  GY[Nv] vertex coordinate
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * nx | double[Nfaces] | x Components of outward normal vector
 * ny | double[Nfaces] | y Components of outward normal vector
 * sJ | double[Nfaces] | the Jacobi of boundary integral
 *
 * @note
 * To get the outward vector, the vertex must be arranged counterclockwise, or the vector
 * is inward and the Jacobi coefficient will be negative.
 *
 */
void Normals(int Nv, double *GX, double *GY, double *nx, double *ny, double *sJ){
    int f, Nfaces = Nv;
    double x1, x2, y1, y2;

    /* outward vector on each faces */
    for(f=0; f<(Nfaces -1); f++){
        x1 = GX[f]; x2 = GX[f+1];
        y1 = GY[f]; y2 = GY[f+1];
        nx[f] =  (y2-y1);
        ny[f] = -(x2-x1);
    }
    /* the last face */
    x1 = GX[f]; x2 = GX[0];
    y1 = GY[f]; y2 = GY[0];
    nx[f] =  (y2-y1);
    ny[f] = -(x2-x1);

    /* normalize the outward vector */
    for(f=0;f<Nfaces;++f){
        sJ[f] = sqrt(nx[f]*nx[f]+ny[f]*ny[f]);
        nx[f] /= sJ[f];
        ny[f] /= sJ[f];
        sJ[f] /= 2.;
    }
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
 * @param [int]      Np
 * @param [double*]  x[Np] x coordinate of nodes
 * @param [double*]  y[Np] y coordinate of nodes
 * @param [double*]  Dr[Np][Np]
 * @param [double*]  Ds[Np][Np]
 *
 * @return
 * return values：
 * name     | type      | description of value
 * -------- |---------- |----------------------
 * drdx     | double*   | the drdx at each node
 * drdy     | double*   | the drdy at each node
 * dsdx     | double*   | the dsdx at each node
 * dsdy     | double*   | the dsdy at each node
 * J        | double*   | the Jacobi at each node
 *
 *
 */
void GeometricFactors(int Np, const double *x, const double *y,
                      double **Dr, double **Ds,
                      double *drdx, double *dsdx,
                      double *drdy, double *dsdy, double *J){

    double *dxdr, *dxds;
    double *dydr, *dyds;
    int n, m;

    dxdr = (double *)calloc(Np, sizeof(double));
    dxds = (double *)calloc(Np, sizeof(double));
    dydr = (double *)calloc(Np, sizeof(double));
    dyds = (double *)calloc(Np, sizeof(double));

    for (n=0;n<Np;n++){
        for (m=0;m<Np;m++){
            dxdr[n] += Dr[n][m]*x[m];
            dxds[n] += Ds[n][m]*x[m];
            dydr[n] += Dr[n][m]*y[m];
            dyds[n] += Ds[n][m]*y[m];
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
