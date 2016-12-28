# Discontinuous Galerkin Ocean Model

# What is it about

Implement the Nodal DGM (Discontinuous Galerkin Method) for various applications, including

* 2d Convection problem
* 2d Shallow water equation

# Features

* Nodal Discontinuous Galerkin method
* Various element shapes, e.g. triangle, quadrilateral
* MPI parallelization

# How to build

* build third party libraries in `ThirdLibrary` folder and install at `3rdParty/install`
    1. parallel-netcdf-1.7.0
    2. ParMetis-3.1
    3. intel MKL library
* set compiler path as environment variable `MPI_DIR`

# Acknowledgement

This project is inspired from [MIDG](http://www.caam.rice.edu/~timwar/RMMC/MIDG.html) program.
Refer to *Nodal discontinuous Galerkin methods - algorithms, analysis, and applications* for more detail.
