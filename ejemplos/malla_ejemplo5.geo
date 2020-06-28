//Placa 2d conductividad de dos capas
gridsize = 0.10;
Point(1) = {0.0, 0.0, 0, gridsize};
Point(2) = {0.8, 0.0, 0, gridsize};
Point(3) = {1.0, 0.0, 0, gridsize};
Point(4) = {1.0, 1.0, 0, gridsize};
Point(5) = {0.8, 1.0, 0, gridsize};
Point(6) = {0.0, 1.0, 0, gridsize};

Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 6};
Line(4) = {6, 1};

Line(5) = {2, 3};
Line(6) = {3, 4};
Line(7) = {4, 5};

Line Loop(9) = {1, 2, 3, 4};
Line Loop(10) = {5, 6, 7, -2};

//Marcas de las condiciones de borde
Physical Line(1) = {1,3,4,5,7};
Physical Line(2) = {6};


Plane Surface(1) = {9};
Plane Surface(2) = {10};

//Marca de los dominios 
Physical Surface(1) = {1};
Physical Surface(2) = {2};


// //Transfinite surface: (le digo que la malla sea estructurada)
// Transfinite Surface {1,2};
