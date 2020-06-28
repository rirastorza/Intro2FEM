//Gmsh ejemplo con CAD de OpenCASCADE.

SetFactory("OpenCASCADE");

Lx = 30;
Ly = 5;
Lz = 5;
r = 1;
R = 10;

Mesh.CharacteristicLengthMax = 2;

Box(1) = {0, 0, 0, Lx, Ly, Lz};
esf1 = newreg;
Sphere(esf1) = {Lx-2, Ly/2, Lz/2, r, -Pi/2, Pi/2, 2*Pi};
esf2 = newreg;
Sphere(esf2) = {Lx-6, Ly/2, Lz/2, r, -Pi/2, Pi/2, 2*Pi};
//Diferencia booleana
BooleanDifference{ Volume{1}; Delete; }{ Volume{esf1}; Volume{esf2}; Delete;}
esf3 = newreg;
Sphere(esf3) = {0, Ly/2, Lz/2, R, -Pi/2, Pi/2, 2*Pi};
//Unión
BooleanUnion{ Volume{1}; Delete; }{ Volume{31}; Delete; }

