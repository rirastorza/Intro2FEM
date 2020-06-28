#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ejemplo: Ecuacion de Poisson 1D pero resuelto con FEniCS 

  -d2u/dx2 = 1   0 < x < L
    u(0) = 0 y u'(L) = constante en los bordes.
    
"""

from __future__ import print_function
from fenics import *
# Defino la malla
nx = 10 #numero de intervalos
l0 = 1.0 #Longitud 1 m
mg = 2.0 #10N
A = 0.01*0.01 #sección 1cmx1cm
E = 210.0e9 #210 GPa
rho = 7850.0 #kg/m**3
vol = A*l0
g = 9.81 #m/s**2
minx, maxx= 0.0, l0 
mesh = IntervalMesh(nx, minx, maxx)#malla en 1D 
V0 = FunctionSpace(mesh, 'CG',1)#Continuous Lagrange elements

# Defino la condición de borde de Dirichlet
def borde_Ar(x, on_boundary): 
    tol = 1.E-14
    return on_boundary and near(x[0], 0.0, tol)


class borde_Ab(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        l0 = 1.0
        return on_boundary and near(x[0], l0, tol)



#Función de la malla con parámetro que indica la topología
marcador_borde = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
#Marcadores en el borde de abajo
b_ab = borde_Ab()
b_ab.mark(marcador_borde, 20)

#Defino el subdominio de los bordes
ds = Measure('ds', domain=mesh, subdomain_data=marcador_borde)
bc_ar = DirichletBC(V0, Constant(0.0), borde_Ar)

bc = [bc_ar]

# Comienzo la formulacion variacional
u = TrialFunction(V0)
v = TestFunction(V0)
f = Constant(rho*vol*g/l0)

#Definicion abstracta 
a = A*E*dot(grad(u), grad(v))*dx#o inner
L = f*v*dx+(mg)*v*ds(20)

# Resuelvo
u = Function(V0)
solve(a == L, u, bc)

print('Tipo de variable:',type(u))

import matplotlib.pyplot as plt

#Extraigo los datos de la solucion u.
uh = u.compute_vertex_values(mesh) 

print('Cantidad de celdas:',nx)
print('Cantidad de vertices:',len(uh))
print('Masa barra: ', rho*vol*g)
print('Solución analítica: ',(mg*l0/(A*E)))

fig, axs = plt.subplots(1,1)

import numpy as np

xu = np.linspace(0.0, 1.0, len(uh),endpoint = True)
xt = np.linspace(0.0, 1.0, 200,endpoint = True)
ut = rho*g*(l0*xt-xt*xt/2.0)/E

axs.plot(xu,uh,'ro',markersize=5)
axs.plot(xt,ut,'b')
axs.hlines((mg*l0/(A*E))+ut[-1],0,l0,linestyles='dashed')
plt.title('Campo de desplazamientos u(x)')


plt.show()

