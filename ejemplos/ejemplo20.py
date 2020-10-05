"""
Ejemplo20.py modificado del tutorial de FEniCS:

https://fenicsproject.org/pub/tutorial/html/._ftut1009.html

Fluido incompresible utilizando Navier-Stokes para un canal (ley de Poisseuille)
en 2D utilizando Incremental Pressure Correction Scheme (IPCS)

  u' + u . nabla(u)) - div(sigma(u, p)) = f
                                 div(u) = 0
"""

from __future__ import print_function
from fenics import *
import numpy as np
import matplotlib.pyplot as plt

T = 10.0           # tiempo final
num_steps = 500    # número de pasos
dt = T / num_steps # paso de tiempo
mu = 1             # viscosidad
rho = 1            # densidad

# Creo la malla
mesh = UnitSquareMesh(16, 16)
V = VectorFunctionSpace(mesh, 'P', 2)
Q = FunctionSpace(mesh, 'P', 1)

#Condiciones de borde
inflow  = 'near(x[0], 0)'#región de presión en la entrada
outflow = 'near(x[0], 1)'#región de presión en la salida
walls   = 'near(x[1], 0) || near(x[1], 1)'#región de velocidad en las paredes

bcu_noslip  = DirichletBC(V, Constant((0, 0)), walls)#velocidad en las paredes
bcp_inflow  = DirichletBC(Q, Constant(8), inflow)#presión en la entrada
bcp_outflow = DirichletBC(Q, Constant(0), outflow)#presión en la salida
bcu = [bcu_noslip]
bcp = [bcp_inflow, bcp_outflow]

# Defino las funciones test y trial
u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

# Defino las funciones test y trial para los pasos de tiempo
u_n = Function(V)
u_  = Function(V)
p_n = Function(Q)
p_  = Function(Q)

# Defino algunas expresiones
U   = 0.5*(u_n + u)
n   = FacetNormal(mesh)
f   = Constant((0, 0))
k   = Constant(dt)
mu  = Constant(mu)
rho = Constant(rho)

# Tensor de deformación
def epsilon(u):
    return sym(nabla_grad(u))

# Tensor de esfuerzo
def sigma(u, p):
    return 2*mu*epsilon(u) - p*Identity(len(u))

# Paso (1) de la formulación variacional
F1 = rho*dot((u - u_n) / k, v)*dx + \
     rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx \
   + inner(sigma(U, p_n), epsilon(v))*dx \
   + dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds \
   - dot(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Paso (2) de la formulación variacional
a2 = dot(nabla_grad(p), nabla_grad(q))*dx
L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx

# Paso (3) de la formulación variacional
a3 = dot(u, v)*dx
L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx

# Junto las matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Aplico las condiciones de borde a las matrices
[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]

# Resuelvo los pasos de tiempo
t = 0
for n in range(num_steps):

    # Actualizo el tiempo
    t += dt

    # Paso (a): velocidad tentativa
    b1 = assemble(L1)
    [bc.apply(b1) for bc in bcu]
    solve(A1, u_.vector(), b1)

    # Paso (b): Corrección de la presión
    b2 = assemble(L2)
    [bc.apply(b2) for bc in bcp]
    solve(A2, p_.vector(), b2)

    # Paso (c): Corrección de la velocidad
    b3 = assemble(L3)
    solve(A3, u_.vector(), b3)

    
    

    # Calculo el error
    u_e = Expression(('4*x[1]*(1.0 - x[1])', '0'), degree=2)
    u_e = interpolate(u_e, V)
    error = np.abs(u_e.vector().get_local() - u_.vector().get_local()).max()
    print('t = %.2f: error = %.3g' % (t, error))
    print('max u:', u_.vector().get_local().max())

    # Actualizo la solución anterior
    u_n.assign(u_)
    p_n.assign(p_)
    
    
plot(u_)
plt.show()


