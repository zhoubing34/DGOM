
Point(1) = {-4000, -4000, 0, 200};
Point(2) = {4000, -4000, 0, 200};
Point(3) = {4000, 4000, 0, 200};
Point(4) = {-4000, 4000, 0, 200};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(3) = {1};

ne = 100;
Transfinite Line{1} = ne+1;
Transfinite Line{2} = ne+1;
Transfinite Line{3} = ne+1;
Transfinite Line{4} = ne+1;

Transfinite Surface{3} = {1,2,4,3};

Physical Line(4) = {3, 4, 1, 2};
Physical Surface(1) = {3};

Mesh.Smoothing = 10;