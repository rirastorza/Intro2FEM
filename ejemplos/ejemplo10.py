#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Ejemplo 5: Ecuacion de Laplace. 

"""

from fenics import *

# Creo la malla y defino la funcion de espacio
mesh = UnitSquareMesh(5, 5)
V = FunctionSpace(mesh, 'CG',1)#Continuous Lagrange elements

# Defino las condiciones de borde
def borde_D(x, on_boundary): 
    tol = 1.E-14
    return on_boundary and near(x[0], 1., tol)

bc_der=DirichletBC(V, Constant(5.0), borde_D)

def borde_I(x, on_boundary):
    tol = 1.E-14
    return on_boundary and near(x[0], 0., tol)

bc_iz = DirichletBC(V, Constant(0.0), borde_I)

def borde_AR(x, on_boundary):
    tol = 1.E-14
    return on_boundary and near(x[1], 0., tol)

bc_arriba = DirichletBC(V, Constant(0.0), borde_AR)

def borde_AB(x, on_boundary):
    tol = 1.E-14
    return on_boundary and near(x[1], 1., tol)

bc_abajo = DirichletBC(V, Constant(0.0), borde_AB)

bc = [bc_iz, bc_der, bc_arriba, bc_abajo]#

#Defino la formulacion variacional
u = TrialFunction(V) 
v = TestFunction(V)
f = Constant(0.0)

#FEniCS detecta solo la forma bilineal y lineal
F = f*v*dx+dot(grad(u), grad(v))*dx
a = lhs(F)
L = rhs(F) 


#####Prueba para obtener matriz de rigidez A
##### https://fenicsproject.org/qa/11091/assemble-a-bilinear-form-returns-a-vector/
A = assemble(a)
A_array = A.array()
print(A_array)

#Resuelvo
u = Function(V)
solve(a == L, u, bc)

vtkfile_T = File('datos/temperatura.pvd')
vtkfile_T << u

#Se puede calcular el flujo de calor (campo vectorial)
V = u.function_space() #Tomo el espacio donde vive u
mesh = V.mesh()#su malla tambien
degree = V.ufl_element().degree() #el elemento finito
W = VectorFunctionSpace(mesh, 'P', degree) #y creo un espacio con los mismos elementos de u

grad_u = project(grad(u), W) #Realizo la proyeccion
k = Constant(1.0)
flux_u = project(-k*grad(u), W)

vtkfile_Q = File('datos/calor.pvd')
vtkfile_Q << flux_u

import matplotlib.pyplot as plt

plt.figure(1)
plot(flux_u)

plt.figure(2)
plot(u.root_node(),mode='warp')

plt.figure(3)
plot(mesh)

plt.show()

#Lista de solvers y precondicionadores 
#list_linear_solver_methods()
#list_krylov_solver_preconditioners()

#Imprime el solver: en este caso PETSc built in LU solver
print(parameters['linear_algebra_backend'])







