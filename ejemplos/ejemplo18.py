#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ejemplo18: Análisis Modal

Extraído y modificado de Jeremy :

https://comet-fenics.readthedocs.io/en/latest/demo/modal_analysis_dynamics/cantilever_modal.py.html

"""
from dolfin import *
import numpy as np

L, B, H = 20., 0.5, 1.

Nx = 200
Ny = int(B/L*Nx)+1
Nz = int(H/L*Nx)+1

mesh = BoxMesh(Point(0.,0.,0.),Point(L,B,H), Nx, Ny, Nz)

#Prueba con Poisson nulo
E, nu = Constant(1e5), Constant(0.)
rho = Constant(1e-3)

#Coeficientes de Lamé
mu = E/2./(1+nu)
lmbda = E*nu/(1+nu)/(1-2*nu)

def eps(v):
    return sym(grad(v))
def sigma(v):
    dim = v.geometric_dimension()
    return 2.0*mu*eps(v) + lmbda*tr(eps(v))*Identity(dim)


V = VectorFunctionSpace(mesh, 'Lagrange', degree=1)
z_ = TrialFunction(V)
v = TestFunction(V)


def izquierda(x, on_boundary):
    return near(x[0],0.)

bc = DirichletBC(V, Constant((0.,0.,0.)), izquierda)

k_form = inner(sigma(v),eps(z_))*dx
l_form = Constant(1.)*z_[0]*dx
K = PETScMatrix()
b = PETScVector()
assemble_system(k_form, l_form, bc, A_tensor=K, b_tensor=b)

m_form = rho*dot(v,z_)*dx
M = PETScMatrix()
assemble(m_form, tensor=M)


eigensolver = SLEPcEigenSolver(K, M)
eigensolver.parameters['problem_type'] = 'gen_hermitian'
eigensolver.parameters['spectral_transform'] = 'shift-and-invert'
eigensolver.parameters['spectral_shift'] = 0.


N_eig = 6
print("Computing {} first eigenvalues...".format(N_eig))
eigensolver.solve(N_eig)

# Exact solution computation
from scipy.optimize import root
from math import cos, cosh

falpha = lambda x: cos(x)*cosh(x)+1
alpha = lambda n: root(falpha, (2*n+1)*pi/2.)['x'][0]

file = File("data/desplazamiento.pvd")

# Extraction
for i in range(N_eig):
    # Extract eigenpair
    r, c, rx, cx = eigensolver.get_eigenpair(i)

    # 3D eigenfrequency
    freq_3D = sqrt(r)/2/pi

    # Beam eigenfrequency
    if i % 2 == 0: # exact solution should correspond to weak axis bending
        I_bend = H*B**3/12.
    else:          #exact solution should correspond to strong axis bending
        I_bend = B*H**3/12.
    rho1 = 1e-3
    E1 = 1e5
    freq_beam = alpha(i/2)**2*sqrt(E1*I_bend/(rho1*B*H*L**4))/2/pi
    
    #print('Solid FE: ',freq_3D,' [Hz]\n')
    #print('Beam theory: ',freq_beam, '[Hz]\n')
    print("Solid FE: {0:8.5f} [Hz]   Beam theory: {1:8.5f} [Hz]".format(freq_3D, freq_beam))

    # Initialize function and assign eigenvector
    eigenmode = Function(V,name="Eigenvector "+str(i))
    eigenmode.vector()[:] = rx
    
    file << eigenmode
    
    
