# Discontinuous Galerkin Ocean Model

![](http://images.cnblogs.com/cnblogs_com/li12242/812136/o_icon_s.png)

# What is it about

Implement the Nodal DGM (Discontinuous Galerkin Method) in various applications, include

* 2d Maxwell equations
* 2d Convection problem
* 2d Shallow water equation
* 3d Ocean modelling (developing)

# Features

* Nodal Discontinuous Galerkin method
* Various element shapes, e.g. triangle, quadrilateral in two-dimensions and prism in three-dimensions (developing)
* Parallel computation

# How to build

* build third party libraries in `ThirdLibrary` folder and install at `3rdParty/install`
    1. parallel-netcdf-1.7.0
    2. ParMetis-3.1
    3. intel MKL library
* set compiler path as environment variable `MPI_DIR`

# Acknowledgement

This project is inspired from [MIDG](http://www.caam.rice.edu/~timwar/RMMC/MIDG.html) program.
Refer to *Nodal discontinuous Galerkin methods - algorithms, analysis, and applications* for more detail.
