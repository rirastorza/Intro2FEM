//Cilindro dos capas
gridsize = 0.5e-3;
r1 = 0.75e-3;
r2 = 0.95e-3;
R = 3.8e-3;


Point(1) = {3.5e-3, 3.5e-3, 0, gridsize};
Point(2) = {-r1+3.5e-3, 3.5e-3, 0, gridsize};
Point(3) = {r1+3.5e-3, 3.5e-3, 0, gridsize};
Point(4) = {-r2+3.5e-3, 3.5e-3, 0, gridsize};
Point(5) = {r2+3.5e-3, 3.5e-3, 0, gridsize};
Point(6) = {-R+3.5e-3, 3.5e-3, 0, gridsize};
Point(7) = {R+3.5e-3, 3.5e-3, 0, gridsize};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 2};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 4};
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 6};

Line Loop(7) = {1,2};
Line Loop(8) = {3,4};
Line Loop(9) = {5,6};

Physical Line(10) = {1,2};
Physical Line(20) = {5,6};

Plane Surface(9) = {8,7};
Plane Surface(10) = {9, 8};

Physical Surface(1) = {10};
Physical Surface(2) = {9};

Mesh.CharacteristicLengthMax = 0.2e-3; 
