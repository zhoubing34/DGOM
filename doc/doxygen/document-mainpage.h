/**

********************************************************************************

@mainpage Discontinuous Galerkin Ocean Model

# Introduction

The DGOM project is a new ocean circulation model with nodal discontinuous Galerkin discretization method. It utilizes the latest high order DG methods and the MPI parallelisation framework to accurately simulate ocean circulation with high efficiency. 

To apply the nodal DG methods to ocean simulations, a fundamental object-oriented library is constructed. The library describes a realistic problem with different level of objects, e.g., standard cell, physical gird and fields. Various choices for standard cells are available, e.g., triangle or quadrilateral for 2-D, prism or hexahedron for 3-D problems. The computational grid are divided with ParMetis library, which maintains the load balanced for all the processes, and are communicated with each other through MPI routines. Boundary conditions are alos very easy to implemented by a user-specific ghost edge functions, and has the generality for different problems.

With the new flexible nodal DG library, multiply applications is solved in the DGOM project, including:

- 2d advection-diffusion equation
- 2d shallow water equations

Each application is designed as an unique solver, which needs special input files. Please refer to page @ref solver for details of each solver's information.

---

# Installing DGOM


---

# Usages of DGOM

With the flexible object-oriented programing method, the fundamental nodal DG library uses multi-level objects to describe a realistic problem. 

Name                    |   Description
-----                   |   ----------
@ref dg_cell            |   standard cell
@ref dg_area            |   geometric information of computational area
@ref dg_phys            |   physical fields

To create a new standard cell, user needs to set the maximum degree N and the cell type. The following code creates three different standard cell with maximum degree N = 3,

    int N = 3;
    dg_cell *point = dg_cell_creat(N, POINT);
    dg_cell *line = dg_cell_creat(N, LINE);
    dg_cell *tri = dg_cell_creat(N, TRIANGLE);

The code `POINT`, `LINE`, and `TRIANGLE` are three enumerators of @ref dg_cell_type which is defined in the header file <@ref dg_cell.h>.

---

@date April 2017
@author li12242, Tianjin University
@email li12242@tju.edu.cn

********************************************************************************
*/