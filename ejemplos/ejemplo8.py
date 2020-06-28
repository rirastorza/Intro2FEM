#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ejemplo 3: Ecuacion de Poisson 1D, hecho a mano 

Deformacion de una barra sometida a tensión axial con resortes en las puntas

  -d2u/dx2 = 1   0 < x < 1
    au'(0) = k0(u(0)-g0) y -au'(1) = k1(u(1)-g1) en los bordes.
    
"""

from __future__ import print_function
import numpy as np #importo numpy y lo denomino np

#Parámetros:
r = 1.e-2 #radio del cable (m)
Area = np.pi*(r**2.) #sección transversal (m^2)
I = 10. #corriente constante (A)
rho = 1.72e-8 #resistividad del cobre (ohm m)
k = 80. #conductividad térmica del cobre (W/m K)
a = Area*k
f = rho*I**2./Area

print(a,f)

#Parámetros de condiciones de borde
k0 = 0.0#coeficiente de transferencia térmica (W/m^2 K)
g0 = 333.15 # 60 °C en Kelvin
k1 = 1e6 #infinito, para simular una condición Dirichlet T = g1 = 20 °C
g1 = 293.15 # 20 °C en Kelvin
#q0 y q1 = 0

#Puntos de x0 a xnx
nx = 85 #numero de intervalos
nodos = nx+1 #cantidad de nodos

uh = np.zeros((nx+1,1))
h = 1./(nx)
Apre = (2./h)*np.eye(nx+1) #(n-1)*(n-1)

rows, cols = np.indices((nx+1,nx+1))
row_vals = np.diag(rows, k=-1)
col_vals = np.diag(cols, k=-1)
z1 = np.zeros((nx+1,nx+1))
z1[row_vals, col_vals]=-1./h

row_vals = np.diag(rows, k=1)
col_vals = np.diag(cols, k=1)
z2 = np.zeros((nx+1,nx+1))
z2[row_vals, col_vals]=-1./h

A = a*(Apre+z1+z2) #Matriz de rigidez

A[0,0] = a/h 
A[nx,nx] = a/h

#Ahora lo nuevo
R = np.zeros((nx+1,nx+1))
R[0,0] = k0
R[nx,nx] = k1


b = f*h*np.ones((nx+1,1))#vector
b[0] = f*h/2.
b[nx] = f*h/2.

#Parte nueva
r = np.zeros((nx+1,1))#vector
r[0] = k0*g0
r[nx] = k1*g1

#Calculo la solución
uh = np.linalg.solve(A+R, b+r)

import matplotlib.pyplot as plt

fig, axs = plt.subplots(1,1)

import numpy as np
#psinew =np.append([0],np.reshape(psi, nx-1))
#uh =np.append(psinew,[0])
#print(psinew)

xu = np.linspace(0, 1.0, nx+1,endpoint = True)


axs.plot(xu,uh,'ro',markersize=5)

print(A)
print(R)
print(b)
print(r)

plt.show()

