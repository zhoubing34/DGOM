/** 

********************************************************************************
@page   overview    Overview

@section Overview

The DGOM project is designed to implement the high order nodal discontinuous Galerkin methods into solving various fluid problems. The nodal DG methods consists of nodal polynomials on standard cell and 

@section Library

With the flexible object-oriented programing method, the fundamental nodal DG library uses multi-level objects to describe a realistic problem. 

Name                    |   Description
-----                   |   ----------
@ref dg_cell            |   standard cell
@ref dg_area            |   geometric information of computational area
@ref dg_phys            |   physical fields

@subsection dg_cell

dg_cell describes the standard cell. It contains three sub-object, some matrixs in user-defined precision and a method to project vertex value to nodes.

Name                    |   Description
-----                   |   ----------
@ref dg_cell_info       |   basic information of standard cell
@ref dg_cell_volume     |   volume information
@ref dg_cell_face       |   face information
f_Dr, f_Ds, f_Dt, f_LIFT|   derivative and lift matrixs
proj_vert2node          |   pronject values from vertexs to nodes

Generate a new dg_cell object, we need to specific the maximum degree N and the cell type. The following code creates three different standard cell with maximum degree N = 3. The code `POINT`, `LINE`, and `TRIANGLE` are three enumerators of @ref dg_cell_type which is defined in the header file <@ref dg_cell.h>.

    int N = 3;
    dg_cell *point = dg_cell_creat(N, POINT);
    dg_cell *line = dg_cell_creat(N, LINE);
    dg_cell *tri = dg_cell_creat(N, TRIANGLE);

@subsection dg_area

dg_area stores the information of geometry grid in parallel. The grid data can be decomposed into four pars, and each part is described by an object,

Name                    |   Description
-----                   |   ----------
@ref dg_grid            |   connection information of each cell
@ref dg_region          |   geometry information of each cell
@ref dg_mesh            |   fetch buffers on mesh with other processes
@ref dg_edge            |   information of global edges

The dg_area is either created by a set of grid files or constructued with uniform grids.


    dg_area *dg_area_create_from_file(dg_cell *cell, char *casename);
    dg_area *dg_area_create_uniform(dg_cell *cell, int Mx, int My,
                                    double xmin, double xmax,
                                    double ymin, double ymax, int type);

For the format of input file, please refer to <@ref dg_grid_reader.c> for details.


@section Solvers

********************************************************************************

@page library:dg_cell library:dg_cell 

The dg_cell class describe the standard element.

To create a new kind of standard cell, user needs to set following basic informations, 

- basic geometry informations
- interpolation nodes
- orthogonal polynomial
- derivative of orthogonal polynomial
- function to project from vertex to nodes

Take triangle as an example, The user needs to specific the basic geometry infomations, e.g., the number of vertexes and faces, the coordinate of vertex and the cell type of each faces. Then need to construct a set of interpolation nodes and the orthogonal polynomials and their derivatives. The project function uses the linear map relation which is the same as coordinate project function.

********************************************************************************

@page library:dg_area library:dg_area

********************************************************************************

@page library:dg_grid library:dg_grid 

The dg_grid class stands for geometry grid. 

It stores the local elements' information and the global vertex coordinates. To create a new dg_grid structure, user can either call the @ref dg_grid_create or the @ref dg_grid_create_from_file2d() function.

The properties of dg_grid includes 

Properties Name     |   Description
-----               |   ----------
Nv                  |   number of vertex on global grid
K                   |   number of element on local grid
EToV                |   local element to vertex list (C style)
EToE                |   adjacent element index (C style)
EToF                |   adjacent face index (C style), e.g., [0, 1, 2] for triangle
EToP                |   process id (C style) of adjacent element 
EToBS               |   face type of each face 
EToR                |   region type of local element 
vx, vy, vz          |   coordinates of the whole vertex

The index of dg_grid, such as the index of elements, vertexs and faces, are stored in the C style, which starts from 0. When these properties were assigned by users directly, please remember to change to 

The methods of dg_grid includes

Methods Name            |   Description
-----                   |   ----------
@ref set_EToBS(dg_grid *grid, int Nsurface, int **SFToV): void |   function to set the EToBS properties
@ref proj_vert2node(dg_grid *grid, int Nfield, double *vertval, double *nodeval): void  |   project the vertex value to the nodes

In the function of set_EToBS, the 

********************************************************************************

@page library:dg_region library:dg_region

********************************************************************************

@page library:dg_mesh library:dg_mesh

********************************************************************************

@page library:dg_edge library:dg_edge

*/