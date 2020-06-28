#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ejemplo: Ecuacion de Poisson 1D pero resuelto con FEniCS 

  -d2u/dx2 = 1   0 < x < 1
    u(0) = 0 y u(1) = 0 en los bordes.
    
"""

from __future__ import print_function
from fenics import *
# Defino la malla
nx = 5 #numero de intervalos
minx, maxx = 0.0, 1.0 
mesh = IntervalMesh(nx, minx, maxx)#malla en 1D 
V = FunctionSpace(mesh, 'P',1)#Lagrange Finite Element

# Defino las condiciones de borde
def borde_D(x, on_boundary): #retorna un boolean
    tol = 1.E-14
    return on_boundary and near(x[0], 1., tol)

def borde_I(x, on_boundary):
    tol = 1.E-14
    return on_boundary and near(x[0], 0., tol)

bc_der = DirichletBC(V, Constant(0.0), borde_D)
bc_iz = DirichletBC(V, Constant(0.0), borde_I)

bc = [bc_iz, bc_der]

# Comienzo la formulacion variacional
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(1.0)
#Definicion abstracta 
a = dot(grad(u), grad(v))*dx #o inner
L = f*v*dx

# Resuelvo
u = Function(V)
solve(a == L, u, bc)

print('Tipo de variable:',type(u))

import matplotlib.pyplot as plt

#Extraigo los datos de la solucion u.
uh = u.compute_vertex_values(mesh) 

print('Cantidad de celdas:',nx)
print('Cantidad de vertices:',len(uh))

fig, axs = plt.subplots(1,1)

import numpy as np

xu = np.linspace(0.0, 1.0, len(uh),endpoint = True)

axs.plot(xu,uh,'ro',markersize=10)

#Comparo con solucion exacta
xe = np.arange(0.0,1.0,0.001)
ue = -0.5*xe*(xe-1.)
axs.plot(xe,ue,'b')

##Tambien se puede calcular en los mismos puntos que uh (para calcular errores)

class Resultado(UserExpression):
    def eval(self, values, x):
        values[0] = -0.5*x[0]*(x[0]-1.0)

u_D = Resultado(degree=1)
u_De = u_D.compute_vertex_values(mesh)

##error_max = np.max(np.abs(u_De - uh))

# Calcula el error en la norma L2
error_L2 = errornorm(u_D, u, 'L2')

##print('Error maximo:',error_max)
print('Error en L2:',error_L2)


axs.plot(xu,u_De,'.b',markersize=10)
plt.title('Soluciones comparadas')


plt.show()

