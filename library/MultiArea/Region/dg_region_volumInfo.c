//
// Created by li12242 on 17/3/7.
//

#include "dg_region_volumInfo.h"

/**
 * @brief calculate the volume factors (drdx, dsdx, drdy and dsdy) for 3d region
 * @param[in,out] region multi-regions object
 * @todo
 */
void dg_reg_volumInfo3d(dg_region *region){

    return;
}

/**
 * @brief calculate the volume factors (drdx, dsdx, drdy and dsdy) for 2d region (triangle and quadrilateral)
 * @details
 * Get the geometric factor in region, the derived geometric factor on each node is calculated as
 * \f[ \left. \frac{\partial x}{\partial r} \right|_{\mathbf{r}_i} =
 * \left. \frac{\partial \mathbf{\varphi}^T_j}{\partial r} \right|_{\mathbf{r}_i} \mathbf{x}_j \f]
 * where \f$ \mathbf{\varphi} \f$ is the vector of Lagrangian basis functions. The same as
 * \f$ \frac{\partial x}{\partial s}, \frac{\partial y}{\partial r} \frac{\partial y}{\partial s} \f$
 *
 * The Jacobian matrix is used to transfer the differentials of (x,y) to (r,s) with the chain rule:
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
 * @param[in,out] region multi-regions object
 */
void dg_reg_volumInfo2d(dg_region *region){

    dg_cell *cell = dg_region_cell(region);
    const int Np = dg_cell_Np(cell);
    double **Dr = dg_cell_Dr(cell);
    double **Ds = dg_cell_Ds(cell);
    const int K = dg_grid_K(region->grid);
    double dxdr[Np], dxds[Np], dydr[Np], dyds[Np];

    // allocation
    region->J = matrix_double_create(K, Np);
    region->drdx = matrix_double_create(K, Np);
    region->drdy = matrix_double_create(K, Np);
    region->dsdx = matrix_double_create(K, Np);
    region->dsdy = matrix_double_create(K, Np);

    double **drdx = region->drdx;
    double **drdy = region->drdy;
    double **dsdx = region->dsdx;
    double **dsdy = region->dsdy;

    int k,n;
    for(k=0;k<K;k++){

        // calculate the dxdr,dxds,dydr,dyds
        matrix_multiply(Np, Np, 1, Dr[0], region->x[k], dxdr);
        matrix_multiply(Np, Np, 1, Ds[0], region->x[k], dxds);
        matrix_multiply(Np, Np, 1, Dr[0], region->y[k], dydr);
        matrix_multiply(Np, Np, 1, Ds[0], region->y[k], dyds);

        for(n=0;n<Np;n++){
            double jtemp = -dxds[n]*dydr[n] + dxdr[n]*dyds[n];

            if(jtemp<0){
                printf("%s (%d): Warning: J[%d][%d] = %lg\n",
                       __FUNCTION__, __LINE__, k, n, region->J[k][n]);
            }
            region->J[k][n] = jtemp;

            /* inverted Jacobian matrix for coordinate mapping */
            drdx[k][n] =  dyds[n]/jtemp; // drdx
            drdy[k][n] = -dxds[n]/jtemp; // drdy
            dsdx[k][n] = -dydr[n]/jtemp; // dsdx
            dsdy[k][n] =  dxdr[n]/jtemp; // dsdy
        }
    }
    return;
}