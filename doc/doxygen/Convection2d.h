/** \page Convection2d
 * Two dimensional scalar convection problem
 *
 * \section Introduction
 * Two dimensional linear transport problem in order to assess
 * the performance of the nodal DG method on triangular elements and on quadrilateral elements.
 *
 * The governing equation is giving as
 * \f[\frac{\partial u}{\partial t} + \nabla \cdot \mathbf{f}(u, x, t) = 0, \quad x\in\Omega \f]
 *
 * where *u* is the scalar variable, the flux term is
 * \f[\mathbf{f} = \mathbf{a}(x) = (a_1(x)u, a_2(x)u).\f]
 *
 * 1. computation domain
 * \f$[-L, L]\f$, where \f$ L = 1 \f$.
 *
 * 2. velocity field
 * \f$[a_1, a_2] = (-wy, wx)\f$, where \f$w = 5\pi/6 \f$. The advection period is \f$T = \frac{2\pi}{w} = 2.4\f$ s.
 * \image html Convect2d_vel.png "The constant velocity field"
 *
 * 3. initial condition
 * \f$ u(\mathbf{x}, t = 0) = exp\left( -\sigma \left| \mathbf{x} - \mathbf{x}_c \right|^2 \right) \f$,
 * where \f$ \mathbf{x}_c = (0, 3/5) \f$, \f$ \sigma = 125 \times 1000/33^2 \f$.
 * \image html Convect2d_initcon.png "The initial scalar distribution"
 *
 * 4. computation meshes
 * different resolutions: [20 x 20], [40 x 40], [60 x 60], [80 x 80].
 * different mesh types:
 *  1. triangle mesh
 *  2. square mesh
 *  3. skewed-rectangular mesh
 *
 * 5. boundary condition
 * transmissive boundary condition (zero gradient).
 *
 * \section Numerical accuracy
 *
 * \section Computation efficiency
 */