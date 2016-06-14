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
 * where *u* is the scalar variable, the flux term is giving by
 * \f[\mathbf{f} = \mathbf{a}(x) = (a_1(x)u, a_2(x)u).\f]
 *
 * Some configuration of the convection problem is given blow.
 *
 * 1. The computation domain is set to \f$[-L, L]\times[-L, L]\f$, where \f$ L = 1 \f$.
 *
 * 2. The const velocity field is \f$[a_1, a_2] = (-wy, wx)\f$ as shown in Fig.1,
 * where \f$w = 5\pi/6 \f$. The period of the field is \f$T = \frac{2\pi}{w} = 2.4\f$ s.
 * \image html Convect2d_vel.png "Fig 1. The constant velocity field"
 *
 * 3. The initial scalar field is giving as \f$ u(\mathbf{x}, t = 0) =
 * exp \left( -\sigma \left| \mathbf{x} - \mathbf{x}_c \right|^2 \right) \f$,
 * where \f$ \mathbf{x}_c = (0, 3/5) \f$, \f$ \sigma = 125 \times 1000/33^2 \f$.
 * \image html Convect2d_initcon.png "Fig 2. The initial scalar distribution"
 *
 * 4. Four different mesh resolutions of uniform meshes are used to compare the performance
 * of different element types, triangle and quadrilateral. The resolution of each mesh is [Ne x Ne],
 * where Ne is the number of elements on each edge and Ne is set to 20, 40, 60 and 80 respectively.
 * The triangle element is derived by subdividing the square element.
 *
 * 5. The Cauchy boundary condition \f$ u(\partial \Omega) =0 \f$ is applied.
 *
 * 6. The time step \f$\Delta t \f$ is obtained by
 * \f[ \Delta t \le CFL \left( \frac{2}{3}min \Delta r_i \right) min_{\Omega}
 * \left( \frac{r_D}{\left| \mathbf{v} \right|} \right)  \f]
 * where \f$r_D\f$ is the radius of circumscribed circle and \f$\Delta r_i\f$
 * is the minimum spatial distance of all nodes in element \f$\Omega\f$.
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
 * To compare the result of different mesh types, several test cases with different
 * orders and mesh resolutions is calculated. The simulation is performed on Mac Mini
 * (2011 Middle) small desktop computer (dual Intel Core i5 64-bit 2.3GHz processors)
 * and each test case is calculated for 5 times to get the averaged computation time.
 *
 * ## Comparison of triangles and quadrilaterals
 *
 * table 1. The numerical error and accuracy of quadrilateral elements
 *
 * N | Grid | dofs | L2 | order | Linf | order | GFLOPS | time |
 * --- | ---  | --- | ---  | ---  | --- | --- | ---- | --- |
 * 1 | 20x20 | 1600 | 1.931838e-02 | 0.00 | 6.998739e-01 | 0.00 | 7.426380 | 0.066945 |
 * 1 | 40x40 | 6400 | 8.865077e-03 | 0.56 | 3.980578e-01 | 0.41 | 6.787010 | 0.573127 |
 * 1 | 60x60 | 14400 | 4.376720e-03 | 0.87 | 2.215936e-01 | 0.72 | 6.882090 | 1.911796 |
 * 1 | 80x80 | 25600 | 2.314170e-03 | 1.11 | 1.283372e-01 | 0.95 | 6.712480 | 4.631606 |
 * Fitted | / | / | / | 0.75 | / | 0.60 | / | / |
 * 2 | 20x20 | 3600 | 4.381359e-03 | 0.00 | 2.265248e-01 | 0.00 | 5.864820 | 0.598873 |
 * 2 | 40x40 | 14400 | 3.292087e-04 | 1.87 | 2.651298e-02 | 1.55 | 6.029960 | 4.741052 |
 * 2 | 60x60 | 32400 | 6.992582e-05 | 1.91 | 4.775822e-03 | 2.11 | 6.061490 | 15.752820 |
 * 2 | 80x80 | 57600 | 2.653846e-05 | 1.68 | 1.674354e-03 | 1.82 | 5.945050 | 38.339500 |
 * Fitted | / | / | / | 1.85 | / | 1.78 | / | / |
 * 3 | 20x20 | 6400 | 5.216514e-04 | 0.00 | 3.211337e-02 | 0.00 | 10.279500 | 1.575742 |
 * 3 | 40x40 | 25600 | 1.738981e-05 | 2.45 | 1.967430e-03 | 2.01 | 10.205500 | 12.530440 |
 * 3 | 60x60 | 57600 | 3.132281e-06 | 2.11 | 4.191399e-04 | 1.91 | 8.271250 | 45.974320 |
 * 3 | 80x80 | 102400 | 9.718241e-07 | 2.03 | 1.370907e-04 | 1.94 | 11.069200 | 94.453960 |
 * Fitted | / | / | / | 2.27 | / | 1.97 | / | / |
 * 4 | 20x20 | 10000 | 6.564403e-05 | 0.00 | 7.407546e-03 | 0.00 | 10.405900 | 5.252896 |
 * 4 | 40x40 | 40000 | 1.143791e-06 | 2.92 | 1.810193e-04 | 2.68 | 11.398600 | 38.888000 |
 * 4 | 60x60 | 90000 | 1.488359e-07 | 2.51 | 1.949072e-05 | 2.75 | 10.872600 | 134.357600 |
 * 4 | 80x80 | 160000 | 5.845488e-08 | 1.62 | 8.404255e-06 | 1.46 | 11.092400 | 317.492400 |
 * Fitted | / | / | / | 2.58 | / | 2.51 | / | / |
 *
 * table 2. The numerical error and accuracy of triangle elements
 *
 * N | Grid | dofs | L2 | order | Linf | order | GFLOPS | time |
 * --- | ---  | --- | ---  | ---  | --- | --- | ---- | --- |
 * 1 | 20x20 | 2400 | 2.197997e-02 | 0.00 | 6.835139e-01 | 0.00 | 5.533340 | 0.154251 |
 * 1 | 40x40 | 9600 | 1.013364e-02 | 0.56 | 3.856273e-01 | 0.41 | 5.575930 | 1.266590 |
 * 1 | 60x60 | 21600 | 5.099316e-03 | 0.85 | 2.140000e-01 | 0.73 | 5.875230 | 3.998704 |
 * 1 | 80x80 | 38400 | 2.778181e-03 | 1.06 | 1.241815e-01 | 0.95 | 5.511430 | 9.925678 |
 * Fitted | / | / | / | 0.73 | / | 0.60 | / | / |
 * 2 | 20x20 | 4800 | 5.870871e-03 | 0.00 | 2.599863e-01 | 0.00 | 7.248070 | 0.711814 |
 * 2 | 40x40 | 19200 | 5.057946e-04 | 1.77 | 3.597915e-02 | 1.43 | 6.920860 | 5.958858 |
 * 2 | 60x60 | 43200 | 1.243913e-04 | 1.73 | 9.185493e-03 | 1.68 | 7.127310 | 19.608140 |
 * 2 | 80x80 | 76800 | 5.471704e-05 | 1.43 | 4.456997e-03 | 1.26 | 7.037860 | 47.465320 |
 * Fitted | / | / | / | 1.70 | / | 1.49 | / | / |
 * 3 | 20x20 | 8000 | 9.543859e-04 | 0.00 | 5.721039e-02 | 0.00 | 5.657270 | 3.679106 |
 * 3 | 40x40 | 32000 | 4.048853e-05 | 2.28 | 4.290283e-03 | 1.87 | 5.882620 | 28.744380 |
 * 3 | 60x60 | 72000 | 7.808962e-06 | 2.03 | 9.234548e-04 | 1.89 | 5.797470 | 95.787140 |
 * 3 | 80x80 | 128000 | 2.431368e-06 | 2.03 | 3.027320e-04 | 1.94 | 5.632380 | 228.853400 |
 * Fitted | / | / | / | 2.16 | / | 1.89 | / | / |
 * 4 | 20x20 | 12000 | 9.432340e-05 | 0.00 | 1.385123e-02 | 0.00 | 5.578670 | 11.590100 |
 * 4 | 40x40 | 48000 | 3.958856e-06 | 2.29 | 4.211664e-04 | 2.52 | 5.687350 | 87.792680 |
 * 4 | 60x60 | 108000 | 7.109361e-07 | 2.12 | 7.760525e-05 | 2.09 | 5.686540 | 305.143800 |
 * 4 | 80x80 | 192000 | 2.157609e-07 | 2.07 | 2.700090e-05 | 1.83 | 5.721570 | 715.588600 |
 * Fitted | / | / | / | 2.20 | / | 2.27 | / | / |
 *
 *
 * \image html Convec2d_errShapeCompare1.png "Fig 3. Comparison of numerical accuracy."
 *
 * \image html Convec2d_timePerDofsShapeCompare1.png "Fig 4. Comparison of computation efficiency."
 *
 * Fig.3 compares the numerical error of different orders and element types on each resolution.
 * It is clear that with the same number of DOFs, results of the quadrilateral elements gives
 * much less numerical error, and it also obtains a higher order of accuracy. To evaluate the
 * performance of computation efficiency, the computation time (without time spent in I/O,
 * initializing and partitioning of the mesh) over different DOFs is
 * illustrated in Fig.4. It is clear that, for the same order of degree, the time spent per
 * DOFs of quadrilaterals is much less than that of triangles. For different order of degrees,
 * as the higher order has less \f$\Delta t\f$ and need much more time step, so the total computation
 * time is multiplied increasing. Finally, through the above exposition, it is obvious that
 * for the same number of DOFs, the quadrilateral meshes give less numerical error and cost
 * less of computation time, which is both better than that of triangle meshes.
 *
 * ## Parallel computation performance: strong scaling
 *
 * table 3. The computation time of different processes with triangle elements [80x80]
 *
 * N | nproc |Grid | dofs | GFLOPS | time/dofs | time |
 * --- | --- | ---  | --- | ---  | ---  | --- |
 * 1 | 1 | 80x80 | 38400 | 5.665302 | 0.000257 | 9.872322 |
 * 1 | 2 | 80x80 | 38400 | 5.450989 | 0.000134 | 5.127562 |
 * 1 | 4 | 80x80 | 38400 | 2.718542 | 0.000134 | 5.139347 |
 * 1 | 6 | 80x80 | 38400 | 1.690168 | 0.000144 | 5.510465 |
 * 1 | 8 | 80x80 | 38400 | 1.266743 | 0.000144 | 5.517339 |
 * 1 | 10 | 80x80 | 38400 | 0.937883 | 0.000155 | 5.959765 |
 * 2 | 1 | 80x80 | 76800 | 7.013565 | 0.000618 | 47.499300 |
 * 2 | 2 | 80x80 | 76800 | 6.458311 | 0.000336 | 25.789613 |
 * 2 | 4 | 80x80 | 76800 | 3.510055 | 0.000309 | 23.724744 |
 * 2 | 6 | 80x80 | 76800 | 2.276470 | 0.000318 | 24.388096 |
 * 2 | 8 | 80x80 | 76800 | 1.732338 | 0.000313 | 24.037503 |
 * 2 | 10 | 80x80 | 76800 | 1.353736 | 0.000320 | 24.606400 |
 * 3 | 1 | 80x80 | 128000 | 5.797042 | 0.001778 | 227.525500 |
 * 3 | 2 | 80x80 | 128000 | 5.423359 | 0.000950 | 121.601625 |
 * 3 | 4 | 80x80 | 128000 | 2.803871 | 0.000919 | 117.597125 |
 * 3 | 6 | 80x80 | 128000 | 1.882649 | 0.000912 | 116.759833 |
 * 3 | 8 | 80x80 | 128000 | 1.410605 | 0.000913 | 116.874469 |
 * 3 | 10 | 80x80 | 128000 | 1.120721 | 0.000919 | 117.684175 |
 *
 * table 4. The computation time of different processes with quadrilateral elements [80x80]
 *
 * N | nproc |Grid | dofs | GFLOPS | time/dofs | time | speedup |
 * --- | --- | ---  | --- | ---  | ---  | --- | --- |
 * 1 | 1 | 80x80 | 25600 | 6.705807 | 0.000181 | 4.632555 | / |
 * 1 | 2 | 80x80 | 25600 | 5.747897 | 0.000106 | 2.710601 | 1.709051 |
 * 2 | 1 | 80x80 | 57600 | 5.913997 | 0.000466 | 26.841640 | / |
 * 2 | 2 | 80x80 | 57600 | 5.506667 | 0.000251 | 14.457620 | 1.856574 |
 * 3 | 1 | 80x80 | 102400 | 10.933375 | 0.000925 | 94.693025 | / |
 * 3 | 2 | 80x80 | 102400 | 10.158062 | 0.000497 | 50.942125 | 1.858835 |
 * 4 | 1 | 80x80 | 160000 | 11.040900 | 0.001986 | 317.797250 | / |
 * 4 | 2 | 80x80 | 160000 | 10.264762 | 0.001068 | 170.886750 | 1.859695 |
 *
 * table 5. The computation time of different processes with triangle elements [80x80]
 *
 * N | nproc |Grid | dofs | GFLOPS | time/dofs | time | speedup |
 * --- | --- | ---  | --- | ---  | ---  | --- | --- |
 * 1 | 1 | 80x80 | 38400 | 5.665302 | 0.000257 | 9.872322 | / |
 * 1 | 2 | 80x80 | 38400 | 5.450989 | 0.000134 | 5.127562 | 1.925344 |
 * 2 | 1 | 80x80 | 76800 | 7.013565 | 0.000618 | 47.499300 | / |
 * 2 | 2 | 80x80 | 76800 | 6.458311 | 0.000336 | 25.789613 | 1.841800 |
 * 3 | 1 | 80x80 | 128000 | 5.797042 | 0.001778 | 227.525500 | / |
 * 3 | 2 | 80x80 | 128000 | 5.423359 | 0.000950 | 121.601625 | 1.871073 |
 * 4 | 1 | 80x80 | 192000 | 5.651868 | 0.003736 | 717.382000 | / |
 * 4 | 2 | 80x80 | 192000 | 5.335824 | 0.001978 | 379.822250 | 1.888731 |
 *
 * Strong scaling indicates the ability to decrease the total run time for a particular problem,
 * scaling across more cores while keeping the overall problem size fixed. In this section,
 * the finest mesh [80x80] is computed with different number of processes. In table 3, different
 * order of triangle meshes is calculated by 1~10 processes. However, for more than 2 processes,
 * the computation time does not decrease and increase a little as each processes need more time
 * to pass information to each other. In table 4 and 5, the computation time coast by different
 * order of triangle and quadrilateral meshes with 1~2 processes is shown. For different orders
 * and element types, computing with 2 processes can reduce almost half of the computation time
 * to double the computation efficiency. In this problem, the ghost edges is adopted in passing
 * the messages between different processes, which is the most efficient method to parallelize
 * the computation algorithm. For ghost elements method, much more information need to be passed,
 * which will bring bad influence in speeding up the computation.
 *
 * ## Slope limiter
 *
 * table 2. The numerical error and accuracy of quadrilateral elements with slope limiter
 *
 * N | Grid | dofs | L2 | order | Linf | order | GFLOPS | time |
 * --- | ---  | --- | ---  | ---  | --- | --- | ---- | --- |
 * 1 | 20x20 | 1600 | 2.268996e-02 | / | 8.145960e-01 | / | 4.283310 | 0.115193 |
 * 1 | 40x40 | 6400 | 1.399883e-02 | 0.35 | 5.869279e-01 | 0.24 | 4.052860 | 0.892397 |
 * 1 | 60x60 | 14400 | 8.170882e-03 | 0.66 | 3.880519e-01 | 0.51 | 4.221700 | 3.283724 |
 * 1 | 80x80 | 25600 | 4.199159e-03 | 1.16 | 2.249757e-01 | 0.95 | 3.250910 | 7.836932 |
 * 2 | 20x20 | 3600 | 4.020016e-03 | / | 2.229774e-01 | / | 3.823500 | 0.847063 |
 * 2 | 40x40 | 14400 | 3.292060e-04 | 1.81 | 2.651292e-02 | 1.54 | 4.608350 | 6.404866 |
 * 2 | 60x60 | 32400 | 6.992121e-05 | 1.91 | 4.775345e-03 | 2.11 | 4.813310 | 20.413860 |
 * 2 | 80x80 | 57600 | 2.653685e-05 | 1.68 | 1.675248e-03 | 1.82 | 4.846390 | 47.292420 |
 * 3 | 20x20 | 6400 | 5.606127e-04 | / | 3.347719e-02 | / | 8.370840 | 1.970642 |
 * 3 | 40x40 | 25600 | 1.738981e-05 | 2.51 | 1.967430e-03 | 2.04 | 8.841920 | 14.669260 |
 * 3 | 60x60 | 57600 | 3.132281e-06 | 2.11 | 4.191399e-04 | 1.91 | 8.478210 | 52.450080 |
 * 3 | 80x80 | 102400 | 9.718241e-07 | 2.03 | 1.370907e-04 | 1.94 | 8.394480 | 121.530200 |
 * 4 | 20x20 | 10000 | 6.564403e-05 | / | 7.407546e-03 | / | 8.071020 | 6.182372 |
 * 4 | 40x40 | 40000 | 1.143791e-06 | 2.92 | 1.810193e-04 | 2.68 | 9.843210 | 45.265460 |
 * 4 | 60x60 | 90000 | 1.488359e-07 | 2.51 | 1.949072e-05 | 2.75 | 9.490330 | 156.578000 |
 * 4 | 80x80 | 160000 | 5.845488e-08 | 1.62 | 8.404255e-06 | 1.46 | 9.315910 | 373.678400 |
 *
 * table 4. The numerical error and accuracy of triangle elements with slope limiter
 *
 * N | Grid | dofs | L2 | order | Linf | order | GFLOPS | time |
 * --- | ---  | --- | ---  | ---  | --- | --- | ---- | --- |
 * 1 | 20x20 | 2400 | 3.106533e-02 | / | 9.239498e-01 | / | 3.190680 | 0.268266 |
 * 1 | 40x40 | 9600 | 2.698678e-02 | 0.10 | 8.547034e-01 | 0.06 | 2.957540 | 1.998792 |
 * 1 | 60x60 | 21600 | 2.414764e-02 | 0.14 | 7.806035e-01 | 0.11 | 3.782030 | 6.605186 |
 * 1 | 80x80 | 38400 | 2.173766e-02 | 0.18 | 7.138205e-01 | 0.16 | 3.595510 | 15.936660 |
 * 2 | 20x20 | 4800 | 2.142572e-02 | / | 9.133729e-01 | / | 4.415130 | 1.188514 |
 * 2 | 40x40 | 19200 | 1.861570e-02 | 0.10 | 8.445063e-01 | 0.06 | 4.320420 | 9.758352 |
 * 2 | 60x60 | 43200 | 2.265943e-04 | 5.44 | 7.378355e-02 | 3.01 | 4.455070 | 29.736800 |
 * 2 | 80x80 | 76800 | 2.148492e-04 | 0.09 | 5.354179e-02 | 0.56 | 5.078560 | 65.839440 |
 * 3 | 20x20 | 8000 | 1.709316e-02 | / | 9.193774e-01 | / | 4.317350 | 4.932582 |
 * 3 | 40x40 | 32000 | 1.098869e-04 | 3.64 | 1.421714e-02 | 3.01 | 4.689970 | 33.784760 |
 * 3 | 60x60 | 72000 | 2.104496e-04 | -0.80 | 4.384325e-02 | -1.39 | 4.732510 | 118.497000 |
 * 3 | 80x80 | 128000 | 4.027016e-04 | -1.13 | 1.061878e-01 | -1.54 | 4.820670 | 273.797800 |
 * 4 | 20x20 | 12000 | 3.479455e-02 | / | 2.213470e+00 | / | 4.417810 | 10.614286 |
 * 4 | 40x40 | 48000 | 1.509942e-04 | 3.92 | 9.642249e-02 | 2.26 | 5.086560 | 100.142980 |
 * 4 | 60x60 | 108000 | 4.506685e-04 | -1.35 | 6.314611e-02 | 0.52 | 4.986920 | 343.553000 |
 * 4 | 80x80 | 192000 | 7.537352e-04 | -0.89 | 1.138729e-01 | -1.02 | 4.996500 | 810.102600 |
 *
 * \image html Convect2d_errLimiterCompare1.png "Fig 5. Comparison of numerical accuracy on different meshes: (a). quadrilateral; (b). triangle."
 *
 * \image html Convect2d_ComEfficientEffect1.png "Fig 6. The effect of slope limiter on computational efficiency of different degrees (quadrilateral elements)"
 *
 *
 */