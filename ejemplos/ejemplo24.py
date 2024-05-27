#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Ejemplo para simular variables/físicas totalmente acopladas

Problema: Ecuación de Poisson mixta.

Última versión:

https://fenicsproject.org/docs/dolfin/2019.1.0/python/demos/mixed-poisson/demo_mixed-poisson.py.html

"""


from dolfin import *
import matplotlib.pyplot as plt


# Create mesh
mesh = UnitSquareMesh(32, 32)

# .. index::
#    pair: FunctionSpace; Brezzi-Douglas-Marini
#    pair: FunctionSpace; Discontinous Lagrange


# Defino los elementos finitos como espacios mixtos 
BDM = FiniteElement("BDM", mesh.ufl_cell(), 1) #Brezzi-Douglas-Marini
DG  = FiniteElement("DG", mesh.ufl_cell(), 0) #Discontinous Lagrange
W = FunctionSpace(mesh, BDM * DG) 


# Defino funciones trial y test
(sigma, u) = TrialFunctions(W)
(tau, v) = TestFunctions(W)

# Defino la función fuente
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)


# Formulación variacional
a = (dot(sigma, tau) + div(tau)*u + div(sigma)*v)*dx
L = - f*v*dx

#Defino la condicion de borde (nunca lo hicimos así)
class BoundarySource(UserExpression):
    def __init__(self, mesh, **kwargs):
        self.mesh = mesh
        super().__init__(**kwargs)
    def eval_cell(self, values, x, ufc_cell):
        cell = Cell(self.mesh, ufc_cell.index)
        n = cell.normal(ufc_cell.local_facet)
        g = sin(5*x[0])
        values[0] = g*n[0]
        values[1] = g*n[1]
    def value_shape(self):
        return (2,)

G = BoundarySource(mesh, degree=2)

# Defino los bordes
def boundary(x):
    return x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS

bc = DirichletBC(W.sub(0), G, boundary) #Constante el flujo!

# Computo la solución
w = Function(W)
solve(a == L, w, bc)
(sigma, u) = w.split() #Separo cada una de las variables

# Plot sigma and u
plt.figure()
plot(sigma)

plt.figure()
plot(u)

plt.show()
