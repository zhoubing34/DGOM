Point(1) = {-300, -300, 0}; Point(2) = {-300, 300, 0};
Point(3) = {300, 300, 0};  Point(4) = {300, -300, 0};
Line(1) = {1, 4}; Line(2) = {4, 3};
Line(3) = {3, 2}; Line(4) = {2, 1};

Line Loop(4) = {1, 2, 3, 4}; 
Plane Surface(10) = {4};
Field[1] = MathEval;
Field[1].F = "50*(x*x + y*y) + 0.5";
Background Field = 1;
Mesh.RecombineAll = 1;
Mesh.Algorithm = 8; // DelQuad (experimental)Coherence;
