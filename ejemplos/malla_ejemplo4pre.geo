//Puntos
Point(1) = {-100, -100, 0, 1};
Point(2) = {100, -300, 0, 1};
Point(3) = {100, 300, 0, 1};
Point(4) = {-100, 100, 0, 1};
Point(5) = {300, 0, 0, 1};
Point(6) = {100, 0, 0, 1};
Point(7) = {-100, 0, 0, 1};
Line(1) = {7, 4};
Line(2) = {4, 3};
Ellipse(3) = {3, 6, 5, 5};
Ellipse(4) = {5, 6, 2, 2};
Line(5) = {2, 1};
Line(6) = {1, 7};
Line(7) = {7, 5};
Line Loop(9) = {2, 3, -7, 1};
Ruled Surface(9) = {9};
Line Loop(11) = {7, 4, 5, 6};
Ruled Surface(11) = {11};
Transfinite Line {1, 3, 6, 4} = 10 Using Progression 1;
Transfinite Line {2, 7, 5} = 10 Using Progression 1;

Transfinite Surface {9};
Transfinite Surface {11};
Recombine Surface {9, 11};
