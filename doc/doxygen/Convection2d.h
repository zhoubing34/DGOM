/** \page Convection2d
 * Two dimensional scalar convection problem
 *
 * \section Introduction
 * A two dimensional linear advection problem is calculated in order to assess
 * the performance of the nodal DG method on triangular and quadrilateral elements.
 *
 * The governing equation is giving as
 * \f[
 * \frac{\partial u}{\partial t} + \nabla \cdot \mathbf{f}(u, x, t) = 0,
 * \quad x\in\Omega
 * \f]
 *
 * Some configuration
 *
 * where *u* is the scalar variable, the flux term is giving by
 * \f[\mathbf{f} = \mathbf{a}(x) = (a_1(x)u, a_2(x)u).\f]
 *
 * 1. computation domain
 * \f$[-L, L]\f$, where \f$ L = 1 \f$.
 *
 * 2. the velocity field
 * \f$[a_1, a_2] = (-wy, wx)\f$, where \f$w = 5\pi/6 \f$. The advection period is
 * \f$T = \frac{2\pi}{w} = 2.4\f$ s.
 * \image html Convect2d_vel.png "Fig 1. The constant velocity field"
 *
 * 3. initial condition
 * The initial condition is giving as \f$ u(\mathbf{x}, t = 0) =
 * exp \left( -\sigma \left| \mathbf{x} - \mathbf{x}_c \right|^2 \right) \f$,
 * where \f$ \mathbf{x}_c = (0, 3/5) \f$, \f$ \sigma = 125 \times 1000/33^2 \f$.
 * \image html Convect2d_initcon.png "Fig 2. The initial scalar distribution"
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
 * \section Discretization
 *
 * To allow information transfer between cells, the Lax-Friedrichs flux is applied
 * in the surface integration as the numerical flux.
 *
 * A slope limiter is adopted in the aim to suppress the numerical oscillation
 * and preserve the numerical stability. However, most slope limiter will introduce
 * too much numerical diffusion
 *
 * \section Result
 *
 * table 1. The original numerical error and accuracy of quadrilateral
 * elements without slope limiter
 *
 * N | Grid | dofs | L2 | order | Linf | order | GFLOPS | time |
 * --- | ---  | --- | ---  | ---  | --- | --- | ---- | --- |
 * 1 | 20x20 | 1600 | 1.931838e-02 | 0.00 | 6.998739e-01 | 0.00 | 7.426380 | 0.066945 |
 * 1 | 40x40 | 6400 | 8.865077e-03 | 0.56 | 3.980578e-01 | 0.41 | 6.787010 | 0.573127 |
 * 1 | 60x60 | 14400 | 4.376720e-03 | 0.87 | 2.215936e-01 | 0.72 | 6.882090 | 1.911796 |
 * 1 | 80x80 | 25600 | 2.314170e-03 | 1.11 | 1.283372e-01 | 0.95 | 6.712480 | 4.631606 |
 * 2 | 20x20 | 3600 | 4.383970e-03 | 0.00 | 2.265248e-01 | 0.00 | 5.471670 | 0.712898 |
 * 2 | 40x40 | 14400 | 3.292089e-04 | 1.87 | 2.651298e-02 | 1.55 | 5.423050 | 5.157016 |
 * 2 | 60x60 | 32400 | 6.992121e-05 | 1.91 | 4.775345e-03 | 2.11 | 5.422430 | 17.694240 |
 * 2 | 80x80 | 57600 | 2.653685e-05 | 1.68 | 1.675248e-03 | 1.82 | 5.271300 | 43.664340 |
 * 3 | 20x20 | 6400 | 5.216514e-04 | 0.00 | 3.211337e-02 | 0.00 | 10.279500 | 1.575742 |
 * 3 | 40x40 | 25600 | 1.738981e-05 | 2.45 | 1.967430e-03 | 2.01 | 10.205500 | 12.530440 |
 * 3 | 60x60 | 57600 | 3.132281e-06 | 2.11 | 4.191399e-04 | 1.91 | 8.271250 | 45.974320 |
 * 3 | 80x80 | 102400 | 9.718241e-07 | 2.03 | 1.370907e-04 | 1.94 | 11.069200 | 94.453960 |
 * 4 | 20x20 | 10000 | 6.564403e-05 | 0.00 | 7.407546e-03 | 0.00 | 10.405900 | 5.252896 |
 * 4 | 40x40 | 40000 | 1.143791e-06 | 2.92 | 1.810193e-04 | 2.68 | 11.398600 | 38.888000 |
 * 4 | 60x60 | 90000 | 1.488359e-07 | 2.51 | 1.949072e-05 | 2.75 | 10.872600 | 134.357600 |
 * 4 | 80x80 | 160000 | 5.845488e-08 | 1.62 | 8.404255e-06 | 1.46 | 11.092400 | 317.492400 |
 *
 * table 2. The numerical error and accuracy of quadrilateral elements with
 * slope limiter
 *
 * N | Grid | dofs | L2 | order | Linf | order | GFLOPS | time |
 * --- | ---  | --- | ---  | ---  | --- | --- | ---- | --- |
 * 1 | 20x20 | 1600 | 2.268996e-02 | 0.00 | 8.145960e-01 | 0.00 | 4.283310 | 0.115193 |
 * 1 | 40x40 | 6400 | 1.399883e-02 | 0.35 | 5.869279e-01 | 0.24 | 4.052860 | 0.892397 |
 * 1 | 60x60 | 14400 | 8.170882e-03 | 0.66 | 3.880519e-01 | 0.51 | 4.221700 | 3.283724 |
 * 1 | 80x80 | 25600 | 4.199159e-03 | 1.16 | 2.249757e-01 | 0.95 | 3.250910 | 7.836932 |
 * 2 | 20x20 | 3600 | 4.020016e-03 | 0.00 | 2.229774e-01 | 0.00 | 3.823500 | 0.847063 |
 * 2 | 40x40 | 14400 | 3.292060e-04 | 1.81 | 2.651292e-02 | 1.54 | 4.608350 | 6.404866 |
 * 2 | 60x60 | 32400 | 6.992121e-05 | 1.91 | 4.775345e-03 | 2.11 | 4.813310 | 20.413860 |
 * 2 | 80x80 | 57600 | 2.653685e-05 | 1.68 | 1.675248e-03 | 1.82 | 4.846390 | 47.292420 |
 * 3 | 20x20 | 6400 | 5.606127e-04 | 0.00 | 3.347719e-02 | 0.00 | 8.370840 | 1.970642 |
 * 3 | 40x40 | 25600 | 1.738981e-05 | 2.51 | 1.967430e-03 | 2.04 | 8.841920 | 14.669260 |
 * 3 | 60x60 | 57600 | 3.132281e-06 | 2.11 | 4.191399e-04 | 1.91 | 8.478210 | 52.450080 |
 * 3 | 80x80 | 102400 | 9.718241e-07 | 2.03 | 1.370907e-04 | 1.94 | 8.394480 | 121.530200 |
 * 4 | 20x20 | 10000 | 6.564403e-05 | 0.00 | 7.407546e-03 | 0.00 | 8.071020 | 6.182372 |
 * 4 | 40x40 | 40000 | 1.143791e-06 | 2.92 | 1.810193e-04 | 2.68 | 9.843210 | 45.265460 |
 * 4 | 60x60 | 90000 | 1.488359e-07 | 2.51 | 1.949072e-05 | 2.75 | 9.490330 | 156.578000 |
 * 4 | 80x80 | 160000 | 5.845488e-08 | 1.62 | 8.404255e-06 | 1.46 | 9.315910 | 373.678400 |
 *
 * table 3. The numerical error and accuracy of triangle elements with
 * slope limiter
 *
 * N | Grid | dofs | L2 | order | Linf | order | GFLOPS | time |
 * --- | ---  | --- | ---  | ---  | --- | --- | ---- | --- |
 * 1 | 20x20 | 2400 | 3.106533e-02 | 0.00 | 9.239498e-01 | 0.00 | 3.190680 | 0.268266 |
 * 1 | 40x40 | 9600 | 2.698678e-02 | 0.10 | 8.547034e-01 | 0.06 | 2.957540 | 1.998792 |
 * 1 | 60x60 | 21600 | 2.414764e-02 | 0.14 | 7.806035e-01 | 0.11 | 3.782030 | 6.605186 |
 * 1 | 80x80 | 38400 | 2.173766e-02 | 0.18 | 7.138205e-01 | 0.16 | 3.595510 | 15.936660 |
 * 2 | 20x20 | 4800 | 2.142572e-02 | 0.00 | 9.133729e-01 | 0.00 | 4.415130 | 1.188514 |
 * 2 | 40x40 | 19200 | 1.861570e-02 | 0.10 | 8.445063e-01 | 0.06 | 4.320420 | 9.758352 |
 * 2 | 60x60 | 43200 | 2.265943e-04 | 5.44 | 7.378355e-02 | 3.01 | 4.455070 | 29.736800 |
 * 2 | 80x80 | 76800 | 2.148492e-04 | 0.09 | 5.354179e-02 | 0.56 | 5.078560 | 65.839440 |
 * 3 | 20x20 | 8000 | 1.709316e-02 | 0.00 | 9.193774e-01 | 0.00 | 4.317350 | 4.932582 |
 * 3 | 40x40 | 32000 | 1.098869e-04 | 3.64 | 1.421714e-02 | 3.01 | 4.689970 | 33.784760 |
 * 3 | 60x60 | 72000 | 2.104496e-04 | -0.80 | 4.384325e-02 | -1.39 | 4.732510 | 118.497000 |
 * 3 | 80x80 | 128000 | 4.027016e-04 | -1.13 | 1.061878e-01 | -1.54 | 4.820670 | 273.797800 |
 * 4 | 20x20 | 12000 | 3.479455e-02 | 0.00 | 2.213470e+00 | 0.00 | 4.417810 | 10.614286 |
 * 4 | 40x40 | 48000 | 1.509942e-04 | 3.92 | 9.642249e-02 | 2.26 | 5.086560 | 100.142980 |
 * 4 | 60x60 | 108000 | 4.506685e-04 | -1.35 | 6.314611e-02 | 0.52 | 4.986920 | 343.553000 |
 * 4 | 80x80 | 192000 | 7.537352e-04 | -0.89 | 1.138729e-01 | -1.02 | 4.996500 | 810.102600 |
 *
 * ### Numerical accuracy
 *
 * \image html Convect2d_errLimiterCompare1.png "Fig 3. The effect of slope limiter on numerical error and accuracy on different meshes (quadrilateral elements)."
 *
 * \image html Convec2d_errShapeCompare1.png "Fig 4. Comparison of the numerical error and accuracy of triangle and quadrilateral elements on different meshes with slope limiter."
 *
 * ### Computation efficiency
 *
 * \image html Convect2d_ComEfficientEffect1.png "Fig 5. The effect of slope limiter on computational efficiency of different degrees (quadrilateral elements)"
 *
 * \image html Convec2d_timePerDofsShapeCompare1.png "Fig 6. Time-per-dof cost of quadrilateral and triangle on different meshes at various orders N."
 *
 * ### Parallel computation performance
 *
 * \image html .png "Fig 7. Total execution time over different number of processes:
 * (a). quadrilateral elements; (b). triangle elements"
 *
 */