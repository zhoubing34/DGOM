[------------]DGOM: 2d swe solver
[------------]case info (1 parameter):
[------------]    1. case name;
ParabolicBowl_Tri
[------------]standard element info (2 parameters):
[------------]    1. element type  |--- 2. triangle
[------------]                     |--- 3. quadrilateral
[------------]    2. order of polynomial
2
1
[------------]mesh info (1 parameters):
[------------]    1. open boundary file;

[------------]time info (2 parameters):
[------------]    1. CFL number;
[------------]    2. dt;
[------------]    3. final time;
0.3
0.0001
1773.1
[------------]physical parameter (3 parameters)
[------------]    1. gravity acceleration;
[------------]    2. hcrit -- minimum water depth of wet cell;
[------------]    3. manning coefficient file;
[------------]    4. initial condition file (h, hu, hv, b);
9.81
1e-2
ParabolicBowl.fric
ParabolicBowl.init
[------------]information for output (2 parameters)
[------------]    1. filename;
[------------]    2. output time interval;
ParabolicBowl_Tri
12
