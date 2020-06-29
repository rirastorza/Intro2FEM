#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ejemplo del libro pag. 37


"""
from __future__ import print_function
from fenics import *
# Defino la malla
nx = 25 #numero de intervalos
minx, maxx = 2.0, 8.0 
mesh = IntervalMesh(nx, minx, maxx)#malla en 1D 
V = FunctionSpace(mesh, 'CG',1)#Continuous Lagrange elements

def borde_iz(x, on_boundary):
    tol = 1.E-14
    return on_boundary and near(x[0], 2., tol)

class borde_der(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        return on_boundary and near(x[0], 8, tol)

bc_iz = DirichletBC(V, Constant(7.0), borde_iz)
bc = [bc_iz]

#Función de la malla con parámetro que indica la topología
marcadores_bordes = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
bcder = borde_der()
bcder.mark(marcadores_bordes, 20)
#Defino el subdominio de los bordes
ds = Measure('ds', domain=mesh, subdomain_data=marcadores_bordes)

# Comienzo la formulacion variacional
u = TrialFunction(V)
v = TestFunction(V)

class funcion_a(UserExpression):
    def eval(self, values, x):
        values[0] = (.5-0.06*x[0])
    def value_shape(self):
        return ()

class funcion_f(UserExpression):
    def eval(self, values, x):
        values[0] = 0.03*pow(x[0]-6,4)
    def value_shape(self):
        return ()
        

a = funcion_a()
f = funcion_f()
g = 0.0#a*p*q

#Definicion abstracta 
a_bilineal = dot(a*grad(u), grad(v))*dx#-a*p*u*v*ds(20)
L_lineal = f*v*dx-g*v*ds(20)

# Resuelvo
u = Function(V)
solve(a_bilineal == L_lineal, u,bc)

print('Tipo de variable:',type(u))

import matplotlib.pyplot as plt

#Extraigo los datos de la solucion u.
uh = u.compute_vertex_values(mesh) 

print('Cantidad de celdas:',nx)
print('Cantidad de vertices:',len(uh))

fig, axs = plt.subplots(1,1)

import numpy as np

xu = np.linspace(2.0, 8.0, len(uh),endpoint = True)

axs.plot(xu,uh,'ro',markersize=5)

plt.show()

