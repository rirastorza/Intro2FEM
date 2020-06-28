//Malla estructurada curva
grilla = 0.1;
Point(1) = {0, 0, 0, grilla};
Point(2) = {-1, -1, 0, grilla};
Point(3) = {-1, 1, 0, grilla};
Point(4) = {-5, -5, 0, grilla};
Point(5) = {-5, 5, 0, grilla};
Line(1) = {4, 2};
Line(2) = {3, 5};
Line(3) = {5, 4};
Circle(4) = {2, 1, 3};
Transfinite Line {3, 4} = 20 Using Progression 1;
Transfinite Line {2, 1} = 40 Using Progression 1;
Line Loop(5) = {3, 1, 4, 2};
Plane Surface(6) = {5};
Transfinite Surface {6};
Recombine Surface{6};


