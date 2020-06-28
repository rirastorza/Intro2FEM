//Ejemplo cilindro 3D de dos capas.
gridsize = .5e-3;
a = 1.5e-3;
b = 1.7e-3;
R = 3.5e-3;
//Circulo del interior
Point(1) = {0, 0, 0, gridsize};
Point(2) = {a, 0, 0, gridsize};
Point(3) = {0, a, 0, gridsize};
Point(4) = {-a, 0, 0, gridsize};
Point(5) = {0,-a, 0, gridsize};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
//Primer capa
Point(6) = {b, 0, 0, gridsize};
Point(7) = {0, b, 0, gridsize};
Point(8) = {-b, 0, 0, gridsize};
Point(9) = {0, -b, 0, gridsize};
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};
//Segunda capa
Point(10) = {R, 0, 0, gridsize};
Point(11) = {0, R, 0, gridsize};
Point(12) = {-R, 0, 0, gridsize};
Point(13) = {0, -R, 0, gridsize};
Circle(9) = {10, 1, 11};
Circle(10) = {11, 1, 12};
Circle(11) = {12, 1, 13};
Circle(12) = {13, 1, 10};
Line Loop(13) = {1,2,3,4};
Line Loop(14) = {5,6,7,8};
Line Loop(15) = {9,10,11,12};
Plane Surface(1) = {13,14};
Plane Surface(2) = {14, 15};
Plane Surface(3) = {13};
Extrude {0, 0, 10.5e-3} {
  Surface{1};
}
Extrude {0, 0, 10.5e-3} {
  Surface{2};
}
Extrude {0, 0, 10.5e-3} {
  Surface{3};
}
Physical Surface(1) = {28,32,36,40};
Physical Surface(2) = {86,90,94,98};
Physical Surface(3) = {44,56};
Physical Volume(1) = {1};
Physical Volume(2) = {2};

