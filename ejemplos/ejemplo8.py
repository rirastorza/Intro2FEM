#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ejemplo: Ecuacion de Poisson 1D pero resuelto con FEniCS 

  -Akd2T/dx2 = f   0 < x < L
    T'(0) = h (T-Tinf) + qinf y T'(L) = h (T-Tinf) + qinf
    
Ejemplo 3.1 de libro: 
https://books.google.com.ar/books/about/Principles_of_Heat_Transfer.html?id=1hVSQBNvr74C&redir_esc=y
"""

from __future__ import print_function
from fenics import *
# Defino la malla
nx = 10 #numero de intervalos

l0 = 0.01 #Longitud 1 cm
A = 0.01*0.1 #sección 1cmx10cm
k = 64.0 # conductividad térmica del acero en [W/mK]
qg = 1.0e6 #En [W/m] 
h = 42.0 #coeficiente de convección en [W/m2K]
Tinf = 353.15 #temperatura en en kelvin

minx, maxx= 0.0, l0 
mesh = IntervalMesh(nx, minx, maxx)#malla en 1D 
V0 = FunctionSpace(mesh, 'CG',1)#Continuous Lagrange elements

class borde_Ar(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        return on_boundary and near(x[0], 0.0, tol)

class borde_Ab(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        l0 = 0.01
        return on_boundary and near(x[0], l0, tol)

#Función de la malla con parámetro que indica la topología
marcador_borde = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
#Marcadores en el borde de abajo
bc_ab = borde_Ab()
bc_ar = borde_Ar()
bc_ab.mark(marcador_borde, 20)
bc_ar.mark(marcador_borde, 30)

#Defino el subdominio de los bordes
ds = Measure('ds', domain=mesh, subdomain_data=marcador_borde)

#bc = [bc_ar,bc_ab]

# Comienzo la formulacion variacional
T = TrialFunction(V0)
v = TestFunction(V0)
f = Constant(qg)

#Definicion abstracta 
a = k*dot(grad(T), grad(v))*dx+h*T*v*ds(20)+h*T*v*ds(30)
L = f*v*dx+h*Tinf*v*ds(20)+h*Tinf*v*ds(30)

# Resuelvo
T = Function(V0)
solve(a == L, T)

import matplotlib.pyplot as plt

#Extraigo los datos de la solucion u.
Th = T.compute_vertex_values(mesh) 

print('Cantidad de celdas:',nx)
print('Cantidad de vertices:',len(Th))

fig, axs = plt.subplots(1,1)

import numpy as np

xu = np.linspace(0.0, 0.01, len(Th),endpoint = True)
axs.plot(xu,Th-273.15,'ro',markersize=5)
axs.set_xlabel('x (m)')
axs.set_ylabel('T (°C)')
plt.title('Temperatura T(x)')


plt.show()

