#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ejemplo16: 

Libro Logg, Mardal, y Wells: Chapter 26. Applications in solid mechanics
Modelo lineal extendido con viscosidad.

"""

from dolfin import *

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True


E = 1000.0
nu = 0.3
mu = E/(2.0*(1.0+nu))
lmbda = E*nu/((1.0+nu)*(1.0-2.0*nu))

#Densidad de masa y coeficiente de viscosidad
rho = 1.0
eta = 500.0

#Parametros de tiempo
beta = 0.25
gamma = 0.5
dt = 0.08
t = 0.0
T = 4.0

#Subdominios
def izquierda(x, on_boundary):
    return near(x[0], 0.) and on_boundary

def derecha(x, on_boundary):
    return near(x[0], 1.) and on_boundary


def actualiza(u,u0,v0,a0,beta,gamma,dt):
    u_vec,u0_vec = u.vector(),u0.vector()
    v0_vec,a0_vec = v0.vector(),a0.vector()
    
    #Actualizar la aceleración y la velocidad
    a_vec = (1.0/(beta*dt*dt))*(u_vec-u0_vec-v0_vec*dt)-(1.0/(2*beta)-1.0)*a0_vec
    
    v_vec = v0_vec+dt*(gamma*a_vec+(1.0-gamma)*a0_vec)
    
    #Acualiza tn <-- tn+1
    v0.vector()[:],a0.vector()[:] = v_vec,a_vec
    u0.vector()[:] = u.vector()


#Malla
mesh = BoxMesh(Point(0., 0., 0.), Point(1., 0.1, 0.04), 60, 10, 5)
V = VectorFunctionSpace(mesh, 'Lagrange', degree=1)

u1 = TrialFunction(V)
w = TestFunction(V)

#Campos del el paso previo (desplazamiento, velocidad, aceleración)
u = Function(V)
u0 = Function(V)
v0 = Function(V)
a0 = Function(V)

boundary_subdomains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundary_subdomains.set_all(0)
force_boundary = AutoSubDomain(derecha)
force_boundary.mark(boundary_subdomains, 3)

#Define la medida de la integral de borde (contorno)
dss = ds(subdomain_data=boundary_subdomains)

#Condición de borde del costado izquierdo
zero = Constant((0.0, 0.0, 0.0))
bc = DirichletBC(V, zero, izquierda)

#Carga aplicada
p0 = 1.
cutoff_Tc = T/5
#NOTAR COMO SE DEFINE LA FUNCIÓN!
h = Expression(("0", "t <= tc ? p0*t/tc : 0","0"), t=0, tc=cutoff_Tc, p0=p0, degree=0)

#Velocidad y aceleración en tn+1
v1=(gamma/(beta*dt))*(u1-u0)-(gamma/beta-1.0)*v0-dt*(gamma/(2.0*beta)-1.0)*a0
a1= (1.0/(beta*dt**2.0))*(u1-u0-dt*v0)-(1.0/(2.0*beta)-1.0)*a0

#Trabajo de fuerzas externas
def Wext(w):
    return dot(w, h)*dss(3)


#Tensor de tensión
def sigma(u,v):
    return 2.0*mu*sym(grad(u))+(lmbda*tr(grad(u))+eta*tr(grad(v)))*Identity(u.geometric_dimension())


#Ecuaciones
F = (rho*dot(a1,w)+inner(sigma(u1,v1),sym(grad(w))))*dx-Wext(w)

a = lhs(F)
L = rhs(F)


#Guardar la solución
file = File("data/desplazamiento.pvd")
u_tip = []
while t <=T:
    t += dt
    h.t = t
    print("Time: ", t)
    #problem.solve()
    solve(a == L, u, bc)
    #u = problem.solve()
    u_tip.append(u(1., 0.05, 0.0)[1])
    actualiza(u,u0,v0,a0,beta,gamma,dt)
    file<<u


# Time-stepping
import matplotlib.pyplot as plt
import numpy as np
#time = np.linspace(0, T, Nsteps+1)
u_final = np.asarray(u_tip)
tiempo = np.arange(0,dt*(len(u_final)),dt)


plt.figure()

plt.plot(tiempo,u_tip)
plt.xlabel("Time")
plt.ylabel("Tip displacement")
#plt.ylim(-0.5, 0.5)
#tiempo2,u_tip2 = np.loadtxt('test.out')
#plt.plot(tiempo2,u_tip2,'-.')

plt.show()
