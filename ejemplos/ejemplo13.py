#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ejemplo13: Pequeña deformación para una barra elástica (voladizo) empotrada en la pared. Ejemplo en 2D.

autor: adaptado de Jeremy Bleyer is a researcher in Solid and Structural Mechanics at Laboratoire Navier, a joint research unit of Ecole Nationale des Ponts et Chaussées, IFSTTAR and CNRS (UMR 8205).

Visitar página:
https://comet-fenics.readthedocs.io/

Resolución teórica de la barra:
http://www.learnaboutstructures.com/Integration-of-the-Curvature-Diagram-to-find-Deflection
"""

from __future__ import print_function
from dolfin import *

L = 25.
H = 1.
Nx = 50
Ny = 10

mesh = RectangleMesh(Point(0., 0.), Point(L, H), Nx, Ny, "crossed")#Malla estructurada

def eps(v):
    return sym(grad(v))

E = Constant(1e5)
nu = Constant(0.3)


mu = E/2/(1+nu)
lmbda = E*nu/(1+nu)/(1-2*nu)
lmbda = 2*mu*lmbda/(lmbda+2*mu)

def sigma(v):
    return lmbda*tr(eps(v))*Identity(2) + 2.0*mu*eps(v)


#Formulación variacional
rho_g = 1e-3
f = Constant((0,-rho_g))

V = VectorFunctionSpace(mesh, 'Lagrange', degree=2)
du = TrialFunction(V)
u_ = TestFunction(V)
a = inner(sigma(du), eps(u_))*dx
l = inner(f, u_)*dx

#Condición de borde
def left(x, on_boundary):
    return near(x[0],0.)

bc = DirichletBC(V, Constant((0.,0.)), left)

#Solución
u = Function(V, name="Desplazamiento")
solve(a == l, u, bc)

plot(1e3*u, mode="Desplazamiento")

# Validación
# ------------------------------
#
# Solution dada por teoría de Euler-Bernoulli.
print("Teoría:", float(3*rho_g*L**4/2/E/H**3))
#Con FEniCS
print("Deflexión máxima:", -u(L,H/2.)[1])


#Posprocesamiento:

Vsig = TensorFunctionSpace(mesh, "DG", degree=0)
sig = Function(Vsig, name="Esfuerzo")
sig.assign(project(sigma(u), Vsig))
print("Esfuerzo en (0,H):", sig(0, H))



import matplotlib.pyplot as plt

plt.figure(1)
plot(1e3*u, mode="displacement")#

plt.figure(2)
plot(mesh)


plt.show()

xdmffile_u = XDMFFile('datos/desplazamiento.xdmf')
xdmffile_u.write(u,0)

## Fields can be exported in a suitable format for vizualisation using Paraview.
## VTK-based extensions (.pvd,.vtu) are not suited for multiple fields and parallel
## writing/reading. Prefered output format is now .xdmf::

#file_results = XDMFFile("resultados_elasticidad.xdmf")
#file_results.parameters["flush_output"] = True
#file_results.parameters["functions_share_mesh"] = True
#file_results.write(u, 0.)
#file_results.write(sig, 0.)
