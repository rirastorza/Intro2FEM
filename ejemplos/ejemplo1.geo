//Ejemplo 1: Placa 3d
Point(1) = {-3.5e-3, 3.5e-3, 0, 1e-3};
Point(2) = {3.5e-3, 3.5e-3, 0, 1e-3};
Point(3) = {3.5e-3, -3.5e-3, 0, 1e-3};
Point(4) = {-3.5e-3, -3.5e-3, 0, 1e-3};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {4, 1, 2, 3};
Plane Surface(6) = {6};
Extrude {0, 0, 1.5e-3} {
 Surface{6};
}
//*Layers{1};
//Recombine;

Physical Volume("interno") = {1};
Physical Surface("frente",1) = {28};
Physical Surface("atras",2) = {6};
Physical Surface("abajo",3) = {27};
Physical Surface("izquierda",4) = {15};
Physical Surface("arriba",5) = {19};
Physical Surface("derecha",6) = {23};
