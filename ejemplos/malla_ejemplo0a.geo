//Placa 2d con cuadrado adentro
gridsize = 0.10;

//Puntos del exterior
Point(1) = {0.0, 0.0, 0, gridsize};
Point(2) = {1.0, 0.0, 0, gridsize};
Point(3) = {1.0, 1.0, 0, gridsize};
Point(4) = {0.0, 1.0, 0, gridsize};

//Puntos del cuadrado interior
Point(5) = {0.25, 0.25, 0, gridsize};
Point(6) = {0.75, 0.25, 0, gridsize};
Point(7) = {0.75, 0.75, 0, gridsize};
Point(8) = {0.25, 0.75, 0, gridsize};

//Lineas exterior
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//Lineas interior 
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line Loop(9) = {1, 2, 3, 4};//Exterior
Line Loop(10) = {5, 6, 7, 8};//Interior

Plane Surface(1) = {9,10};//Superficie con agujero
Plane Surface(2) = {10};

Physical Line(1) = {1,3,4};//borde a temperatura baja
Physical Line(2) = {2};//borde a temperatura alta
// 
Physical Surface(10) = {1};
Physical Surface(20) = {2};
