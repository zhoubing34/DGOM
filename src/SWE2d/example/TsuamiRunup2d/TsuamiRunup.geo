lland = .1;
lopen = .05;
Point(1) = {0, 0, 0, lland};
Point(2) = {5.488, 0, 0, lopen};
Point(3) = {0, 3.402, 0, lland};
Point(4) = {5.488, 3.402, 0, lopen};
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Line(2) = {1, 2, 3};
Physical Line(4) = {4};
Physical Surface(1) = {1};