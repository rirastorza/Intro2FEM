// No estructurada adaptativa
L = 10;
grilla = 1;
grillaFina = .1;
Point(1) = {-L/2, -L/2, 0, grilla};
Point(2) = {L/2, -L/2, 0, grilla};
Point(3) = {L/2, L/2, 0, grillaFina};
Point(4) = {-L/2, L/2, 0, grillaFina};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4}; 	
Plane Surface(6) = {5}; 
