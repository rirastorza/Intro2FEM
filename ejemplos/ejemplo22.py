#!/usr/bin/env python
# coding: utf-8

"""
ejemplo13: 

Libro: Larson y Bengzon, capítulo 11.

y en:

https://comet-fenics.readthedocs.io

"""

from __future__ import print_function
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt

L, H = 5, 0.3
mesh = RectangleMesh(Point(0., 0.), Point(L, H), 100, 10, "crossed")

def laterales(x, on_boundary):
    return (near(x[0], 0) or near(x[0], L)) and on_boundary
def abajo(x, on_boundary):
    return near(x[1], 0) and on_boundary
def arriba(x, on_boundary):
    return near(x[1], H) and on_boundary


VT = FunctionSpace(mesh, "CG", 1)#Lagrange de orden 1 (temperatura es un escalar)
T_, dT = TestFunction(VT), TrialFunction(VT) #Ojo T_ es la función de test (v)
Delta_T = Function(VT, name="Incremento de la temperatura")
aT = dot(grad(dT), grad(T_))*dx
LT = Constant(0)*T_*dx

bcT = [DirichletBC(VT, Constant(50.), abajo), 
       DirichletBC(VT, Constant(0.), arriba),
       DirichletBC(VT, Constant(0.), laterales)]
solve(aT == LT, Delta_T, bcT)#Pone en Delta_T el valor despejado
plt.figure()
p = plot(Delta_T,title="Delta T", mode="contour")
plt.colorbar(p)
#plt.show()


E = Constant(50e3)
nu = Constant(0.2)
mu = E/2/(1+nu)
lmbda = E*nu/(1+nu)/(1-2*nu)
alpha = Constant(1e-5)

rho_g = 1e-3
b = Constant((0, -rho_g))

rho_g = 0.
b = Constant((0, -rho_g))

def eps(v):
    return 0.5*(nabla_grad(v) + nabla_grad(v).T)

def sigma(v, dT):
    return (lmbda*tr(eps(v))- alpha*(3*lmbda+2*mu)*dT)*Identity(2) + 2.0*mu*eps(v)

Vu = VectorFunctionSpace(mesh, 'CG', 2) #Es una tensor
du = TrialFunction(Vu)
u_ = TestFunction(Vu)#Ojo esta es la función test
Wint = inner(sigma(du, Delta_T), eps(u_))*dx #Delta_T ya es  conocida del anterior
aM = lhs(Wint)
LM = rhs(Wint) + inner(b, u_)*dx

bcu = DirichletBC(Vu, Constant((0., 0.)), laterales)

u = Function(Vu, name="Desplazamiento")


solve(aM == LM, u, bcu)

plt.figure()
p = plot(1e3*u,title="Desplazamiento [mm]", mode="displacement")
plt.colorbar(p)
#plt.show()
plt.figure()
p = plot(sigma(u, Delta_T)[0, 0],title="Tensión horizontal [MPa]")
plt.colorbar(p)
plt.show()

vtkfile_u = File('datos/desplazamiento.pvd')
vtkfile_u << u


##Si se agrega el peso
#rho_g = 2400*9.81e-6
#f.assign(Constant((0., -rho_g))) 
#solve(aM == LM, u, bcu)

#plt.figure()
#p = plot(1e3*u[1],title="Desplazamiento vertical [mm]")
#plt.colorbar(p)
#plt.show()
#plt.figure()
#p = plot(sigma(u, Delta_T)[0, 0],title="Tensión horizontal [MPa]")
#plt.colorbar(p)
#plt.show()

