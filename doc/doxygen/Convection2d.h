/** \page Convection2d
 * 二维对流扩散方程求解器
 *
 * \section 简介
 * 二维对流扩散方程作为最典型的守恒方程形式，经常用于检验数值格式的误差，精度以及守恒性等。
 *
 * 守恒形式控制方程为
 * \f[
 * \frac{\partial T}{\partial t} + \frac{\partial E}{\partial x} + \frac{\partial G}{\partial y} =
 * \frac{\partial }{\partial x}\left( \mu_x \frac{\partial T}{\partial x} \right)
 * + \frac{\partial }{\partial y}\left( \mu_y \frac{\partial T}{\partial y} \right),
 * \quad (x,y)\in\Omega
 * \f]
 *
 * 其中 *T* 为守恒变量，而 \f$ \mathbf{F} = \left[ E,G \right] = \left[ uT, vT \right] \f$ 为输运通量项，
 * \f$ \mu=\left[\mu_x, \mu_y \right] \f$ 分别为 *x*,*y* 方向的扩散系数。
 *
 * 验证求解器计算的基本算例包括以下几种。
 *
 * \subsection 旋转流
 *
 * 算例基本设置如下
 *
 * 1. 计算域大小为 \f$[-1, 1]\times[-1, 1]\f$，采用均匀网格计算，可以指定 *x*,*y* 维度上单元个数。
 *
 * 2. 流速场逆时针旋转流动 \f$ \left( u,v \right) = (-wy, wx)\f$ （如图1所示），其中 \f$w = 5\pi/6 \f$，
 * 运动周期为 \f$t = 2\pi/w = 2.4\f$ s。
 * ![](Convection2d/rotation_flowfield.png "图1. 旋转流速场")
 *
 * 3. 初始场采用高斯函数表达形式 \f$ T(\mathbf{x}, t=0) =
 * exp \left[ -\sigma \left| \mathbf{x} - \mathbf{x}_c \right|^2 \right] \f$，其中初始位置
 * \f$ \mathbf{x}_c = (0, 3/5) \f$，参数 \f$ \sigma = 125 \times 1000/33^2 \f$.
 * ![](Convection2d/rotation_initialfield.png "图2. 初始场")
 *
 * 4. 旋转流动算例只是纯对流运动，不包括扩散运动。
 *
 * 5. 开边界条件采用柯西边界条件 \f$ u(\partial \Omega) =0 \f$。
 *
 * 6. 时间步长 \f$\Delta t \f$ 根据如下公式确定
 * \f[ \Delta t = CFL \cdot min_{\Omega_i} \left( \frac{r_D}{| \mathbf{v}| (N+1)} \right) \f]
 * 其中 \f$r_D\f$ 为单元长度，二维单元中采用公式 \f$ r_D = \sqrt{A/\pi} \f$ 计算。*N* 为单元内数值解
 * 的最高阶数，\f$ \mathbf{v} \f$ 为单元内速度矢量大小。
 *
 * \subsection 对流扩散运动
 *
 * 对流扩散运动采用与旋转流基本相同的设置，如计算域、网格等等。其不同设置包括
 *
 * 1. 流速场为常速度场 \f$ u=v=0.5 \f$。
 * 2. 扩散系数全场为常数 \f$ \mu_x = \mu_y = \mu_0 \f$。
 * 3. 初始解与扩散系数有关，其表达式为
 * \f[ T(\mathbf{x}, t=0)= exp \left[ -\frac{(x-x_0)^2}{\mu_0}-\frac{(y-y_0)^2}{\mu_0}\right] \f]
 * 其中初始位置 \f$ (x_0, y_0)=(0.5,0.5) \f$。计算终止时，解析解表达式为
 * \f[ T(\mathbf{x}, t=t_0) = \frac{1}{4t+1} exp \left[ -\frac{(x-x_0-ut_0)^2}{\mu_0(4t+1)}
 * -\frac{(y-y_0-vt_0)^2}{\mu_0(4t+1)} \right] \f]
 *
 * \section 数值离散
 *
 * 数值离散采用两次分部积分的强形式方法计算，具体离散过程如下。
 *
 * 由于方程具有二阶偏导项，首先利用辅助变量 \f$ (p,q) \f$ 转换为一阶偏导形式
 *
 * \f[ \left\{ \begin{array}{l}
 * \frac{\partial T}{\partial t} + \frac{\partial E}{\partial x} + \frac{\partial G}{\partial y} =
 * \frac{\partial p}{\partial x} + \frac{\partial q}{\partial y} \cr
 * p = \mu_0 \frac{\partial T}{\partial x} \cr
 * q = \mu_0 \frac{\partial T}{\partial y} \cr
 * \end{array} \right. \f]
 *
 * 将方程组乘以试验函数后在计算域内积分，并且采用两次分部积分后得到
 *
 * \f[ \left\{ \begin{array}{r}
 * \int_{\Omega_i} l_i \left( \frac{\partial T}{\partial t} + \frac{\partial E}{\partial x}
 * + \frac{\partial G}{\partial y} \right) \mathrm{d\Omega} + \oint_{\partial \Omega_i}
 *  l_i \left( E^* n_x + G^* n_y- En_x - Gn_y \right) \mathrm{\partial \Omega} \cr
 * = \int_{\Omega_i} l_i \left( \frac{\partial p}{\partial x} +
 * \frac{\partial q}{\partial y} \right) \mathrm{d\Omega} + \oint_{\partial \Omega_i}
 * l_i ( p^* n_x+q^* n_y-p n_x-q n_y ) \mathrm{\partial \Omega} \cr
 * \int_{\Omega_i} l_i p \mathrm{d\Omega} = \int_{\Omega_i} l_i \mu_0 \frac{\partial T}{\partial x}
 * \mathrm{d\Omega} + \oint_{\partial \Omega_i} \mu_0 l_i (T^* n_x - T n_x) \mathrm{\partial \Omega} \cr
 * \int_{\Omega_i} l_i q \mathrm{d\Omega} = \int_{\Omega_i} l_i \mu_0 \frac{\partial T}{\partial y}
 * \mathrm{d\Omega} + \oint_{\partial \Omega_i} \mu_0 l_i (T^* n_y - T n_y) \mathrm{\partial \Omega} \cr
 * \end{array} \right. \f]
 *
 * 其中 \f$ E^*, G^* \f$ 为对流项数值通量，\f$ p^*, q^*, T^* \f$ 为扩散项数值通量。
 *
 * 写为矩阵形式
 *
 * \f[ \left\{ \begin{array}{r}
 * J \left( M \frac{\partial T_j}{\partial t} + S_x E_j + S_y G_j \right)
 * + J_s \cdot M_{es} \left( E^*_j n_x + G^*_j n_y - E_j n_x - G_j n_y \right) \cr
 * = J \left( S_x p_j + S_y q_j \right) + J_s \cdot M_{es}
 * \left( p^*_j n_x + q^*_j n_y - p_j n_x - q_j n_y \right) \cr
 * J \cdot M p_j = \mu J \cdot S_x T_j + \mu J_s \cdot M_{es} \left( T^*_j n_x - T_j n_x \right) \cr
 * J \cdot M q_j = \mu J \cdot S_y T_j + \mu J_s \cdot M_{es} \left( T^*_j n_y - T_j n_y \right) \cr
 * \end{array}\right. \f]
 *
 * 由于对流通量中包含流速场，在程序中，将流速场作为守恒变量 *T* 的分量，对应的输运通量分量则全部为0，即
 * \f[ T=\begin{bmatrix} c \cr u \cr v \end{bmatrix},
 * E=\begin{bmatrix} uc \cr 0 \cr 0 \end{bmatrix},
 * G=\begin{bmatrix} vc \cr 0 \cr 0 \end{bmatrix}\f]
 *
 * 输运项数值通量采用迎风通量，其通量形式为
 *
 * \f[ \mathbf{F}^* = \{\{ \mathbf{F} \}\} + \frac{C}{2}[[T]] \f]
 *
 * 其中，\f$ \{\{\mathbf{F}\}\} = 0.5\left( \mathbf{F}^- + \mathbf{F}^+ \right) \f$，
 * \f$ [[T]] = T^- \mathbf{n}^- + T^+ \mathbf{n}^+ \f$，
 * \f$ C=\left| \mathbf{n}\cdot\mathbf{u} \right| \f$。
 *
 * 注意，表达式 \f$ \{\{\mathbf{F}\}\} \f$ 与 \f$ [[T]] \f$ 始终满足数值通量的交换性，
 * 即 \f$ \mathbf{F}^*(U^-,U^+,n^-,n^+)=\mathbf{F}^*(U^+,U^-,n^+,n^-) \f$。
 *
 * 计算粘性项的数值通量则稍微复杂一点，若将 \f$ (p^*, q^*) \f$ 用向量形式 \f$ \mathbf{q}^* \f$ 表示，
 * 那么根据前人研究，\f$ \mathbf{q}^* \f$ 与 \f$ T^* \f$ 表达式有以下几种 (Li, 2005)
 *
 * | Method |     q* |     T* |
 * | :----  | :----: | :----: |
 * | LDG (Cockburn and Shu, 1998) | \f$ \{\{q\}\} - C11[[T]] - C12[[q]] \f$ | \f$ \{\{T\}\}+C12[[T]] \f$ |
 * | DG  (Castillo and Cockburn, 2000) | \f$ \{\{q\}\} - C11[[T]] - C12[[q]] \f$ | \f$ \{\{T\}\}+C12[[T]]-C22[[q]] \f$ |
 * | Brezzi et al. (2000) | \f$ \{\{q\}\} - a[[T]] \f$ | \f$ \{\{T\}\} \f$ |
 * | Bassi-Rebay (1997)   | \f$ \{\{q\}\} \f$ | \f$ \{\{T\}\} \f$ |
 *
 * 其中DG格式中 C12=0.5 时等价于迎风格式，C11 与 C22 为非负实数，用来保证格式稳定性。
 *
 * \section 计算结果
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
 * \section 参考文献
 *
 * 1. Cockburn B, Shu CW. The Local Discontinuous Galerkin Method for
 * Time-Dependent Convection-Diffusion Systems. SIAM J. Numer. Anal. 1998; 35: 2440–2463.
 * 2. Castillo P, Cockburn B, Perugia I, Schotzau D. An a priori Error Analysis of the
 * Local Discontinuous Galerkin Method for Elliptic Problems. SIAM J. Numer. Anal. 2000; 38: 1676–1706.
 * 3. Brezzi F, Manzini G, Marini D, Pietra P, Russo A. Discontinuous Galerkin Approximations
 * for Elliptic Problems. Numer. Methods for Partial Differ. Equ. 2000; 16: 365–378.
 * 4. Bassi F, Rebay S. A High-Order Accurate Discontinuous Finite Element Method for the Numerical
 * Solution of the Compressible Navier–Stokes Equations. J. Comput. Phys. 1997; 131: 267–279.
 * 5. Li B Q. Discontinuous finite elements in fluid dynamics and heat transfer[M].
 * Springer Science & Business Media, 2005.
 */