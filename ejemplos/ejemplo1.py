#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ejemplo 1: polinomios lineales a trozo 

Jupyter notebook
    
autor: Ramiro Irastorza
"""

from __future__ import print_function
import numpy as np #importo numpy y lo denomino np
import matplotlib.pyplot as plt

#Puntos de x0 a xnx
nx = 5 #numero de intervalos
nodos = nx+1 #cantidad de nodos
x = np.linspace(0,1,200) #continuo
v = 2.0*x*np.sin(2.0*np.pi*x)+3.0 #funcion
xi = np.linspace(0,1,nodos) #nodos equiespaciados
vi = 2.0*xi*np.sin(2.0*np.pi*xi)+3.0 #valores en los nodos

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig, axs = plt.subplots(1,1)
axs.plot(xi,vi,'-ro',markersize=10)
axs.plot(x,v,'b')
axs.set_ylim(-1,6)
axs.axhline(0, color='gray')
axs.vlines(xi[0],vi[0],0,linestyles='dashed')
axs.annotate(r'$x_{0}$', xy=(xi[0]-0.02, -0.5),fontsize=16)
axs.vlines(xi[-1],vi[-1],0,linestyles='dashed')
axs.annotate(r'$x_{n}$', xy=(xi[-1]-0.02, -0.5),fontsize=16)

plt.show()

