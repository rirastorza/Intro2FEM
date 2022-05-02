// Malla ejemplo 2D
// the square
Point(1) = {0, 0, 0,0.1};
Point(2) = {0, 1, 0,0.1};
Point(3) = {1, 1, 0,0.1};
Point(4) = {1, 0, 0,0.1};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
// creating the surfaces
Line Loop(6) = {2, 3, 4, 1};
Plane Surface(8) = {6};

Physical Surface(2) = {8};

Physical Line(20) = {4};
