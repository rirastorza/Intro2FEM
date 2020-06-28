//Geometria
Merge "hueso.stl";
//Cilindro
des = 0.1;
radioin = 0.55;
Point(1) = {0+des,0+des,1.6,1};
Point(2) = {2.5+des,0+des,1.6,1};
Point(3) = {0+des,2.5+des,1.6,1};
Point(4) = {-2.5+des,0+des,1.6,1};
Point(5) = {0+des,-2.5+des,1.6,1};
Point(6) = {radioin+des,0+des,1.6,1};
Point(7) = {0+des,radioin+des,1.6,1};
Point(8) = {-radioin+des,0+des,1.6,1};
Point(9) = {0+des,-radioin+des,1.6,1};
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Circle(5) = {6,1,7};
Circle(6) = {7,1,8};
Circle(7) = {8,1,9};
Circle(8) = {9,1,6};
Line Loop(5) = {1,2,3,4};
Line Loop(6) = {5,6,7,8};
Plane Surface(6) = {5};
Plane Surface(7) = {6};
Extrude {0,0,-3.5} {
   Surface{6};
} 
Extrude {0,0,-1.6} {
   Surface{7};
}
Delete{Volume{1}; Surface{6};}
Delete{Volume{2}; Surface{7};}
Line Loop(53) = {5,6,7,8};
Line Loop(54) = {1,2,3,4};
Plane Surface(55) = {53,54};
Surface Loop(2) = {-1,30,21,25,17,29,55,-43,-52,-39,-51,-47};
Volume(3) = {2};
Surface Loop(3) = {1};
Volume(4) = {3};
Physical Surface(20) = {30,25,29,21,17};//Superficie exterior
Physical Surface(30) = {52,47,39,51,43};//Superficie del termistor
Physical Volume(300) = {3};//Médula ósea
Physical Volume(400) = {4};//Matriz ósea
