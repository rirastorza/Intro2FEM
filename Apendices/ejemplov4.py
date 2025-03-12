# -*- coding: utf-8 -*-

#https://github.com/FEniCS/dolfinx/blob/v0.9.0.post1/python/test/unit/fem/test_interpolation.py#L1080
#Lo pregunté acá: https://fenicsproject.discourse.group/t/interpolation-to-and-from-sub-meshes/17182
from mpi4py import MPI

import numpy as np
# import pytest

import basix
import ufl
from basix.ufl import blocked_element, custom_element, element, enriched_element, mixed_element
from dolfinx import default_real_type, default_scalar_type
from dolfinx.fem import (
    Expression,
    Function,
    assemble_scalar,
    create_interpolation_data,
    form,
    functionspace,
)
from dolfinx.geometry import bb_tree, compute_collisions_points
from dolfinx.mesh import (
    CellType,
    create_mesh,
    create_rectangle,
    create_submesh,
    create_unit_cube,
    create_unit_square,
    locate_entities,
    locate_entities_boundary,
    meshtags,
)

# parametrize_cell_types = pytest.mark.parametrize(
#     "cell_type",
#     [
#         CellType.interval,
#         CellType.triangle,
#         CellType.tetrahedron,
#         CellType.quadrilateral,
#         CellType.hexahedron,
#     ],
# )

"""Test interpolation of a function between a sub-mesh and its parent mesh."""
mesh = create_unit_square(MPI.COMM_WORLD, 6, 7)

def left_locator(x):
    return x[0] <= 0.5 + 1e-14

def ref_func(x):
    return x[0] + x[1] ** 2

tdim = mesh.topology.dim
cells = locate_entities(mesh, tdim, left_locator)
submesh, parent_cells, _, _ = create_submesh(mesh, tdim, cells)

u0 = Function(functionspace(mesh, ("Lagrange", 2)))
u0.interpolate(ref_func)

V1 = functionspace(submesh, ("DG", 0))
u1 = Function(V1)

# Interpolate u0 (defined on 'full' mesh) into u0 (defined on
# 'sub'0mesh)
u1.interpolate(u0, cells0=parent_cells, cells1=np.arange(len(parent_cells)))

u1_exact = Function(V1)
u1_exact.interpolate(ref_func)
atol = 5 * np.finfo(default_scalar_type).resolution
# np.testing.assert_allclose(u1_exact.x.array, u1.x.array, atol=atol)

# Map from sub to parent
W = functionspace(mesh, ("DG", 0))
w = Function(W)

# Interpolate Function defined on sub-mesh (u1_exact) to the part of
# a Function on the full mesh (w)
w.interpolate(u1, cells0=np.arange(len(parent_cells)), cells1=parent_cells)
w_exact = Function(W)
w_exact.interpolate(ref_func, cells0=cells)
# np.testing.assert_allclose(w.x.array, w_exact.x.array, atol=atol)


import dolfinx
with dolfinx.io.XDMFFile(mesh.comm, "f_mesh.xdmf", "w") as file:
    file.write_mesh(mesh)
    file.write_function(w)

with dolfinx.io.XDMFFile(submesh.comm, "f_submesh.xdmf", "w") as file:
    file.write_mesh(submesh)
    file.write_function(u1)
