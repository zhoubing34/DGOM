lland = 500;
lopen = 500;
Point(1) = {0, -250, 0, lopen};
Point(5) = {10000, -250, 0, lland};
Point(2) = {20000, -250, 0, lopen};
Point(3) = {20000, 250, 0, lopen};
Point(6) = {10000, 250, 0, lland};
Point(4) = {0, 250, 0, lopen};
Line(1) = {1, 2};
Line(2) = {2, 3}; // east boundary
Line(3) = {3, 4}; // west boundary
Line(4) = {4, 1};
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Line(2) = {1, 3};
Physical Line(5) = {2};
Physical Line(4) = {4};
Physical Surface(1) = {1};