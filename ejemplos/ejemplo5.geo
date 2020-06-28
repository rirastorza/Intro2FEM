//Importo la malla en .stl
Merge "cavidad.stl";
RefineMesh;
CreateTopology;
ll[] = Line "*";
For j In {0 : #ll[]-1}
  Compound Line(newl) = ll[j];
EndFor
ss[] = Surface "*";
s = news;
For i In {0 : #ss[]-1}
  Compound Surface(s+i) = ss[i];
  Printf("Superficie con numero %g",ss[i]) ;
EndFor
Surface Loop(1) = {1};//Cavidad
Mesh.RemeshAlgorithm = 1; 
Mesh.RemeshParametrization = 7;
Geometry.HideCompounds = 0; 
Volume(1) = {1};
