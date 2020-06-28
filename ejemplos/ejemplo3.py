#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ejemplo 3: aproximación de una función

Jupyter notebook
    
autor: Ramiro Irastorza
"""


from __future__ import print_function
import numpy as np #importo numpy y lo denomino np
import matplotlib.pyplot as plt

#Definimos la función sombrero
def phi(x,xi,i):
    hi = xi[1]-xi[0]
    f = np.zeros_like(x)
    for m in range(len(x)):
        if i == 0:
            if xi[i] <= x[m] < xi[i+1]:
                f[m] = (xi[i+1]-x[m])/(xi[i+1]-xi[i])
            else:
                f[m] = 0.0        
        elif i == len(xi):
            if (xi[i-1] < x[m] <= xi[i]):
                f[m] = (x[m]-xi[i-1])/(xi[i]-xi[i-1])
            else:
                f[m] = 0.0
        else:
            if (xi[i-1] < x[m] <= xi[i]):
                f[m] = (x[m]-xi[i-1])/(xi[i]-xi[i-1])
            elif xi[i] < x[m] <= xi[i+1]:
                f[m] = (xi[i+1]-x[m])/(xi[i+1]-xi[i])
            else:
                f[m] = 0.0
        
    return f


#Puntos de x0 a xnx
nx = 5 #numero de intervalos
nodos = nx+1 #cantidad de nodos
x = np.linspace(0,1,200) #continuo
y = 2.0*x*np.sin(2.0*np.pi*x)+3.0 #funcion
xi = np.linspace(0,1,nodos) #nodos equiespaciados
vi = 2.0*xi*np.sin(2.0*np.pi*xi)+3.0 #valores en los nodos

v = vi[0]*phi(x,xi,0)+vi[1]*phi(x,xi,1)+vi[2]*phi(x,xi,2)+vi[3]*phi(x,xi,3)+vi[4]*phi(x,xi,4)+vi[5]*phi(x,xi,5)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig, axs = plt.subplots(1,1)
axs.plot(xi,vi,'ro',markersize=10)
axs.plot(x,v,'r')
axs.plot(x,y,'b')
axs.fill_between(x, y, v)
axs.set_ylim(-1,5)
axs.axhline(0, color='gray')
axs.vlines(xi[0],vi[0],0,linestyles='dashed')
axs.annotate(r'$x_{0}$', xy=(xi[0]-0.02, -0.5),fontsize=16)
axs.vlines(xi[-1],vi[-1],0,linestyles='dashed')
axs.annotate(r'$x_{n}$', xy=(xi[-1]-0.02, -0.5),fontsize=16)

plt.show()



