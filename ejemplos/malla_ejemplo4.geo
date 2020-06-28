// Inputs
L = 2;
grilla = 1;

// Geometry
Point(1) = {-L/2, -L/2, 0, grilla};
Point(2) = {L/2, -L/2, 0, grilla};
Point(3) = {L/2, L/2, 0, grilla};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};

Line Loop(5) = {1, 2, 3}; 	
Plane Surface(6) = {5}; 

// 
Transfinite Line {1,2,3} = 12 Using Progression 1;
Transfinite Surface {6};
// Recombine Surface {6};
