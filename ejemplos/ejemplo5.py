#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ejemplo: Ecuacion de Poisson 1D, hecho a mano 

  -d2u/dx2 = 1   0 < x < 1
    u(0) = 0 y u(1) = 0 en los bordes.
    
"""

from __future__ import print_function
import numpy as np #importo numpy y lo denomino np

#Puntos de x0 a xnx
nx = 5 #numero de intervalos
nodos = nx+1 #cantidad de nodos

uh = np.zeros((nx+1,1))
h = 1./(nx)
Apre =(2./h)*np.eye(nx-1) #(n-1)*(n-1)
#


rows, cols = np.indices((nx-1,nx-1))
row_vals = np.diag(rows, k=-1)
col_vals = np.diag(cols, k=-1)
z1 = np.zeros((nx-1,nx-1))
z1[row_vals, col_vals]=-1./h



row_vals = np.diag(rows, k=1)
col_vals = np.diag(cols, k=1)
z2 = np.zeros((nx-1,nx-1))
z2[row_vals, col_vals]=-1./h



A = Apre+z1+z2 #Matriz de rigidez

print(A)

f = 1.
b = f*h*np.ones((nx-1,1))#vector

print(b)


#Calculo la solución
#A xi = b

ji = np.linalg.solve(A, b)



import matplotlib.pyplot as plt

fig, axs = plt.subplots(1,1)

print(ji)

print(np.reshape(ji, nx-1))

##import numpy as np
jinew =np.append([0],np.reshape(ji, nx-1))
uh =np.append(jinew,[0])
#print(jinew)

xu = np.linspace(0, 1.0, nx+1,endpoint = True)


axs.plot(xu,uh,'rs',markersize=10)

#Comparo con solucion exacta
xe = np.arange(0.0,1.0,0.001)
ue = -0.5*xe*(xe-1.)
axs.plot(xe,ue,'b')

plt.show()

