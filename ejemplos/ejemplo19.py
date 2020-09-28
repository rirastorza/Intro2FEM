#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ejemplo19: Análisis Modal 2D con deformaciones planas


"""
from dolfin import *
import numpy as np
from matplotlib import pyplot as plt

mesh = RectangleMesh(Point(0., 0.), Point(10, 2), 500, 100)
E = 70e9 # modulus [Pa]
nu = 0.33
rho = 2700 # density [kg/m^3]

mu = E/2./(1.+nu) # Lame constant
lmbda = E*nu/((1.+nu)*(1.-2.*nu)) # Lame constant

def eps(v):
    return sym(grad(v))

def sigma(v):
    dim = v.geometric_dimension()
    return 2.0*mu*eps(v) + lmbda*tr(eps(v))*Identity(dim)

V = VectorFunctionSpace(mesh, 'Lagrange', degree=2)
u_ = TrialFunction(V)
du = TestFunction(V)

def left(x, on_boundary):
    return near(x[0], 0.)
bc = DirichletBC(V, Constant((0., 0.)), left)

k_form = inner(sigma(du), eps(u_))*dx
l_form = Constant(1.)*u_[0]*dx
K = PETScMatrix()
b = PETScVector()
assemble_system(k_form, l_form, bc, A_tensor=K, b_tensor=b)

m_form = rho*dot(du, u_)*dx
M = PETScMatrix()
assemble(m_form, tensor=M)
bc.zero(M)#Try zero-ing out the rows of the mass matrix M corresponding to the constrained degrees of freedom,

eigensolver = SLEPcEigenSolver(K, M)
eigensolver.parameters['problem_type'] = 'gen_hermitian'
eigensolver.parameters['spectral_transform'] = 'shift-and-invert'
eigensolver.parameters['spectral_shift'] = 0.

N_eig = 10
print("Computing {} first eigenvalues...".format(N_eig))
eigensolver.solve(N_eig)

for i in range(N_eig):
    r, c, rx, cx = eigensolver.get_eigenpair(i)
    freq_2d = sqrt(r)/2/pi
    print(freq_2d)
    print(c)

    eigenmode = Function(V, name="Eigenvector"+str(i))
    eigenmode.vector()[:] = rx
    plot(eigenmode)
    plt.show()
