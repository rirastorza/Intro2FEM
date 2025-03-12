# -*- coding: utf-8 -*-

# https://fenicsproject.discourse.group/t/problem-transferring-facet-tags-to-submesh/11213/2
# https://fenicsproject.discourse.group/t/extend-submesh-function-to-full-mesh/8006/5


import dolfinx
from dolfinx.io.gmshio import read_from_msh
from mpi4py import MPI
import numpy as np



mesh, cell_tags, facet_tags = dolfinx.io.gmshio.read_from_msh("modelo1.msh", MPI.COMM_WORLD, 0)
#
# print(type(facet_tags))
#
# # Check physical volumes and surfaces of the mesh
# print(np.unique(cell_tags.values))
# print(np.unique(facet_tags.values))

# Create submesh which is the union of physical volumes 200, 300 and 400
submesh, entity_map, vertex_map, geom_map = dolfinx.mesh.create_submesh(mesh, mesh.topology.dim, cell_tags.indices[(cell_tags.values==200)|(cell_tags.values==300)|(cell_tags.values==400)])

# Transfer facet tags from parent mesh to submesh
tdim = mesh.topology.dim
fdim = tdim - 1

c_to_f = mesh.topology.connectivity(tdim, fdim)
f_map = mesh.topology.index_map(fdim)
all_facets = f_map.size_local + f_map.num_ghosts
print(all_facets)
all_values = np.zeros(all_facets, dtype=np.int32)
all_values[facet_tags.indices] = facet_tags.values
print(np.unique(all_values))

#NUEVO
c_to_c = mesh.topology.connectivity(tdim, tdim)
fDOMAIN_map = mesh.topology.index_map(tdim)
all_cells = fDOMAIN_map.size_local + fDOMAIN_map.num_ghosts
print(all_cells)
all_cells = np.zeros(all_cells, dtype=np.int32)
all_cells[cell_tags.indices] = cell_tags.values
print(np.unique(all_cells))

#-----


submesh.topology.create_entities(fdim)
subf_map = submesh.topology.index_map(fdim)
submesh.topology.create_connectivity(tdim, fdim)
c_to_f_sub = submesh.topology.connectivity(tdim, fdim)
num_sub_facets = subf_map.size_local + subf_map.num_ghosts
sub_values = np.empty(num_sub_facets, dtype=np.int32)
for i, entity in enumerate(entity_map):
    parent_facets = c_to_f.links(entity)
    child_facets = c_to_f_sub.links(i)
    for child, parent in zip(child_facets, parent_facets):
        sub_values[child] = all_values[parent]

sub_meshtag = dolfinx.mesh.meshtags(submesh, submesh.topology.dim-1, np.arange(num_sub_facets, dtype=np.int32), sub_values)
submesh.topology.create_connectivity(submesh.topology.dim-1, submesh.topology.dim)
print(type(sub_meshtag))
with dolfinx.io.XDMFFile(MPI.COMM_WORLD, "submesh_facet.xdmf", "w") as xdmf:
    xdmf.write_mesh(submesh)
    xdmf.write_meshtags(sub_meshtag,submesh.geometry)
    # xdmf.write_meshtags(sub_meshtag)

#NUEVO
submesh.topology.create_entities(tdim)
subf_map_cells = submesh.topology.index_map(tdim)
submesh.topology.create_connectivity(tdim, tdim)
c_to_c_sub = submesh.topology.connectivity(tdim, tdim)
num_sub_cells = subf_map_cells.size_local + subf_map_cells.num_ghosts
sub_values_cells = np.empty(num_sub_cells, dtype=np.int32)
for i, entity in enumerate(entity_map):
    parent_cells = c_to_c.links(entity)
    child_cells = c_to_c_sub.links(i)
    for child, parent in zip(child_cells, parent_cells):
        sub_values_cells[child] = all_cells[parent]
#
sub_meshtag_cells = dolfinx.mesh.meshtags(submesh, submesh.topology.dim, np.arange(num_sub_cells, dtype=np.int32), sub_values_cells)
submesh.topology.create_connectivity(submesh.topology.dim, submesh.topology.dim)
# #-----------------------------------------




with dolfinx.io.XDMFFile(MPI.COMM_WORLD, "submesh_cells.xdmf", "w") as xdmf:
    xdmf.write_mesh(submesh)
    xdmf.write_meshtags(sub_meshtag_cells,submesh.geometry)
# print(np.unique(sub_meshtag.values))


print(np.unique(sub_meshtag_cells.values))
print(np.unique(sub_meshtag.values))

#Resolución del problema eléctrico en la submalla.

from dolfinx.fem import functionspace
V = functionspace(submesh, ("Lagrange", 1))

from dolfinx import fem
import numpy as np
from dolfinx import default_scalar_type

bcu_electrode_activo = fem.dirichletbc(default_scalar_type(40),  fem.locate_dofs_topological(V, fdim, sub_meshtag.find(40)), V)
bcu_electrode_pasivo = fem.dirichletbc(default_scalar_type(0),  fem.locate_dofs_topological(V, fdim, sub_meshtag.find(30)), V)

#https://jsdokken.com/dolfinx-tutorial/chapter3/subdomains.html
Q = functionspace(submesh, ("DG", 0))

conductivity = fem.Function(Q)
# ventricle = sub_meshtag_cells.find(100)
# conductivity.x.array[ventricle] = np.full_like(ventricle, 2.0, dtype=default_scalar_type)
electrodes = sub_meshtag_cells.find(200)
conductivity.x.array[electrodes] = np.full_like(electrodes, 4.6e6, dtype=default_scalar_type)
plastic = sub_meshtag_cells.find(300)
conductivity.x.array[plastic] = np.full_like(plastic, 1e-5, dtype=default_scalar_type)
gm = sub_meshtag_cells.find(400)
conductivity.x.array[gm] = np.full_like(gm, 0.152, dtype=default_scalar_type)

import ufl

u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
f = fem.Constant(submesh, default_scalar_type(0))

a = ufl.inner(conductivity*ufl.grad(u), ufl.grad(v)) * ufl.dx
L = f * v * ufl.dx

from dolfinx.fem.petsc import LinearProblem

problem = LinearProblem(a, L, bcs=[bcu_electrode_activo,bcu_electrode_pasivo], petsc_options={"ksp_type": "cg", "pc_type": "hypre","ksp_rtol ": 1e-10 , "ksp_atol ": 1e-15})

uh = problem.solve()

from dolfinx import io
from pathlib import Path
results_folder = Path("results")
results_folder.mkdir(exist_ok=True, parents=True)
filename = results_folder / "tension"

#-----------------
#Solo para dibujar
with io.XDMFFile(submesh.comm, filename.with_suffix(".xdmf"), "w") as xdmf:
    xdmf.write_mesh(submesh)
    xdmf.write_function(uh)

print(type(uh))


#Interpolate in the parent mesh

# Map from sub to parent
W = functionspace(mesh, ("DG", 0))
w = fem.Function(W)

# Interpolate Function defined on sub-mesh (u1_exact) to the part of
# a Function on the full mesh (w)
w.interpolate(uh, cells0=np.arange(len(entity_map)), cells1=entity_map)


import dolfinx
with dolfinx.io.XDMFFile(mesh.comm, "f_parent_mesh.xdmf", "w") as file:
    file.write_mesh(mesh)
    file.write_function(w)
