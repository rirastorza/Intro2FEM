#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Ejemplo 11: FEniCS geometría 2D: cilindros concéntricos
Problema térmico no lineal

"""
from __future__ import print_function

import os

string = "gmsh -2 mallaejemplo11.geo"
os.system(string)

string = "dolfin-convert mallaejemplo11.msh mallaejemplo11.xml"
os.system(string)

from fenics import *

mesh = Mesh("mallaejemplo11.xml");
subdomains = MeshFunction('size_t',mesh,"mallaejemplo11_physical_region.xml");
boundary_markers = MeshFunction('size_t',mesh,"mallaejemplo11_facet_region.xml");

tol = 1E-14
#Constantes termicas
k_0 = Constant(10.0)
k_1 = Constant(400.0)

#Fuente, Temperatura variable
Tb = 310.15 #En Kelvin
DeltaT = 4.0

V = FunctionSpace(mesh, 'CG', 2)

#Defino condiciones de contorno por medio del archivo xml
bx0 = DirichletBC(V, Constant(Tb+DeltaT), boundary_markers,10)
bx1 = DirichletBC(V, Constant(Tb), boundary_markers, 20)
bcs = [bx0,bx1]

class K(UserExpression):
    def __init__(self, subdomains, k_0, k_1, **kwargs):
        super().__init__(**kwargs)
        self.subdomains = subdomains
        self.k_0 = k_0
        self.k_1 = k_1        
    def eval_cell(self, values, x, cell):
        if self.subdomains[cell.index] == 2:
            values[0] = self.k_0
        else:
            values[0] = self.k_1

kappa = K(subdomains, k_0, k_1, degree=2)

## Define variational problem
u = Function(V)  # Note: not TrialFunction!
v = TestFunction(V)
f = Constant(0.00)

dx = dx(subdomain_data=subdomains)
ds = ds(subdomain_data=boundary_markers)
# Define variational problem
F = kappa*dot(grad(u), grad(v))*dx-f*v*dx

solve(F == 0, u, bcs)

#Cálculo del flujo
ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)
n = FacetNormal(mesh)
flux = -(kappa)*dot(grad(u),n)*ds(20)
total_flux = assemble(flux)
print('Flujo numerico:',total_flux)

import numpy as np

#vtkfile_u = File('salidas/temperatura.pvd')
#vtkfile_u << u

#vtkfile_dominios = File('salidas/dominios.pvd')
#vtkfile_dominios << subdomains

#import matplotlib.pyplot as plt

#Comparación con teoría
kb = 400.0
ka = 10.0
r1 = 0.75e-3;
r2 = 0.95e-3;
r3 = 3.8e-3;

Req = ln(r2/r1)/(2.0*np.pi*ka) +ln(r3/r2)/(2.0*np.pi*kb) #Resistencia por unidad de long.

flujoTeorico = DeltaT/Req #flujo por unidad de longitud
print('Flujo teorico:',flujoTeorico)
