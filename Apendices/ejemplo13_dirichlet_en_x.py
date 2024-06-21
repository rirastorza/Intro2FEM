#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ejemplo13: Pequeña deformación en barra de LxH dejando libre uy, pero fijando ux = 0 cuando x = L. Se agregó Dirichlet en la base, es decir, ux = 0 y uy = 0 cuando y = 0.

autor: adaptado de:

https://jsdokken.com/dolfinx-tutorial/chapter3/component_bc.html
"""

from __future__ import print_function
from dolfin import *
L = 1
H = 1.3
lambda_ = 1.25
mu = 1
rho = 1
g = 1
Nx = 30
Ny = 30



mesh = RectangleMesh(Point(0., 0.), Point(L, H), Nx, Ny, "crossed")#Malla estructurada
V = VectorFunctionSpace(mesh, 'Lagrange', degree=1)

def eps(v):
    #return sym(grad(v))
    return 0.5*(nabla_grad(v) + nabla_grad(v).T)

def sigma(v):
    return lambda_*tr(eps(v))*Identity(2) + 2.0*mu*eps(v)


#Condiciones de borde
def abajo(x, on_boundary):
    return near(x[1],0.)

def derecha(x, on_boundary):
    return near(x[0],L)

bc = DirichletBC(V, Constant((0.0,0.0)), abajo)
bcx = DirichletBC(V.sub(0), Constant(0.0), derecha)

bcs = [bc, bcx]

#Formulación variacional
f = Constant((0.0,-rho*g))
T = Constant((0.0,0.0))

du = TrialFunction(V)
u_ = TestFunction(V)

a = inner(sigma(du), eps(u_))*dx
l = inner(f, u_)*dx + inner(T, u_) *ds

#Solución
u = Function(V, name="Desplazamiento")
solve(a == l, u, bcs)

import matplotlib.pyplot as plt

plt.figure(1)
plot(u, mode="displacement")#

plt.show()
