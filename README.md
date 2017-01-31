# Discontinuous Galerkin Ocean Model

# What is it about

Implement the Nodal DGM (Discontinuous Galerkin Method) for various applications, including

* 2d Convection problem
* 2d Shallow water equation

# Features

* Nodal Discontinuous Galerkin method
* Various element shapes, e.g. triangle, quadrilateral
* MPI parallization

# How to build


* build third party libraries in `externlib` folder and install at `external/install`
    1. parallel-netcdf-1.7.0
    2. parMetis-3.1
    3. metis-5.1.0
    4. clapack-3.1.1.1
* set compiler path as environment variable `MPI_DIR`

# Acknowledgement

This project is inspired from [MIDG](http://www.caam.rice.edu/~timwar/RMMC/MIDG.html) program.
Refer to *Nodal discontinuous Galerkin methods - algorithms, analysis, and applications* for more detail.
