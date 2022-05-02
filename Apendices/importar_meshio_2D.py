#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Ejemplo para importar geometria 2D a FEniCS

https://fenicsproject.discourse.group/t/different-results-for-the-same-mesh-created-differently-hyperelasticity/3944/2

Ver también:

https://jsdokken.com/converted_files/tutorial_pygmsh.html#second

"""

from dolfin import *
import meshio
import numpy as np

#Primera parte convierte .msh a xdmf
nombre_malla = 'geometria2D'
msh = meshio.read(nombre_malla+'.msh')

print(msh)
#Levanto los triangulos
points = msh.points
points = points[:,:2]
cells = np.vstack([cells.data for cells in msh.cells if cells.type == "triangle"])
cell_data = np.vstack([msh.cell_data_dict["gmsh:physical"][key] for key in msh.cell_data_dict["gmsh:physical"].keys() if key =="triangle"])
mesh_new = meshio.Mesh(points=points, cells=[("triangle", cells)], cell_data={"name_to_read": cell_data})
meshio.write(nombre_malla+".xdmf", mesh_new)

#Aquí los bordes
facet_cells = np.vstack([cells.data for cells in msh.cells if cells.type == "line"])
facet_data = np.vstack([msh.cell_data_dict["gmsh:physical"][key] for key in msh.cell_data_dict["gmsh:physical"].keys() if key == "line"])
facet_mesh = meshio.Mesh(points=points,cells=[("line", facet_cells)],cell_data={"name_to_read": facet_data})
meshio.write(nombre_malla+"_bordes.xdmf", facet_mesh)


#Ahora pasamos a fenics
mesh = Mesh()

#Ahora levanta el xdmf de los dominios
with XDMFFile(nombre_malla+".xdmf") as infile:
    infile.read(mesh)
    mvcv = MeshValueCollection("size_t", mesh, mesh.topology().dim())
    infile.read(mvcv, "name_to_read")

subdomains = cpp.mesh.MeshFunctionSizet(mesh, mvcv)

mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim()-1)

#Ahora levanta el xdmf de los bordes
with XDMFFile(nombre_malla+"_bordes.xdmf") as infile:
    infile.read(mvc, "name_to_read")

boundaries = cpp.mesh.MeshFunctionSizet(mesh, mvc)

n = FacetNormal(mesh)
V = VectorFunctionSpace(mesh, 'CG', degree=2)


import matplotlib.pyplot as plt

plt.figure(1)
plot(mesh)
plt.figure(2)
plot(subdomains)

plt.show()


