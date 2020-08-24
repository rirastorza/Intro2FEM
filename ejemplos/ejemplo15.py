#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ejemplo 15: Problema dependiente del tiempo ejemplo Difusión del calor. 
Euler hacia adelante.
Solución diseñada exacta.
    
"""

from __future__ import print_function
from fenics import *
import numpy as np

T = 10.0            # tiempo final
Nt = 40000     # número de pasos
dt = T / Nt # paso de tiempo

nx = 20 #numero de intervalos
minx, maxx = 0.0, 1.0 
mesh = IntervalMesh(nx, minx, maxx)#malla en 1D 
V = FunctionSpace(mesh, 'P',1)#Lagrange Finite Element

## Defino las condiciones de borde
beta = 1.0
u_D = Expression('1+x[0]*x[0]+beta*t', degree=1, beta=beta, t=0)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)#Condición de borde
u_n = interpolate(u_D, V)#Condición inicial

## formulación variacional
u = TrialFunction(V)
v = TestFunction(V)
f_n = Constant(beta-2) #Al ser constante, vale lo mismo para todas las muestras

F = u*v*dx + dt*dot(grad(u_n), grad(v))*dx - (u_n + dt*f_n)*v*dx
a, L = lhs(F), rhs(F)


u = Function(V)
u_inter = []
t = 0.0

for nn in range(Nt):
    
          
    t += dt
    u_D.t = t
    
    # calcula la solución
    solve(a == L, u, bc)
    # Temperatura a la mitad de la barra
    u_inter.append(u(0.5))    
    
    #Calcula errores
    u_exacta = interpolate(u_D,V)
    error = np.abs(u_exacta.vector()-u.vector()).max()
    print('t= %.2f: error = %.3g' %(t,error))
    
    
    u_n.assign(u)
    
    #plot(u)
    

print('Tipo de variable:',type(u))

import matplotlib.pyplot as plt
import numpy as np

#Extraigo los datos de la solucion u.
uh = u.compute_vertex_values(mesh) 

print('Cantidad de celdas:',nx)
print('Cantidad de vertices:',len(uh))

xu = np.linspace(0.0, 1.0, len(uh),endpoint = True)

plt.subplot(2,1,1)
plt.plot(xu,uh,'ro',markersize=10)

##Comparo con solucion exacta
xe = np.arange(0.0,1.0,0.001)
ue = 1.0+xe**2.0+beta*T
plt.plot(xe,ue,'b')
plt.ylabel('Solución en t = 10 s')
plt.xlabel('x')

plt.subplot(2,1,2)
tiempo = np.arange(0.0,T,dt)
plt.plot(tiempo,u_inter,'b')
plt.ylabel('Solución en x = 0.5')
plt.xlabel('t (s)')

plt.show()

print(dt)
