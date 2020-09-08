#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ejemplo17: Pequeña deformación dinámica para una barra elástica (voladizo) empotrada en la pared. Ejemplo en 3D.

autor: adaptado de Jeremy Bleyer is a researcher in Solid and Structural Mechanics at Laboratoire Navier, a joint research unit of Ecole Nationale des Ponts et Chaussées, IFSTTAR and CNRS (UMR 8205).

Visitar página:
https://comet-fenics.readthedocs.io/
"""

from dolfin import *
import numpy as np
import matplotlib.pyplot as plt


parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

# Propiedades del material
E  = 1000.0
nu = 0.3
mu    = Constant(E / (2.0*(1.0 + nu)))
lmbda = Constant(E*nu / ((1.0 + nu)*(1.0 - 2.0*nu)))

#------------------------
#Discretización en tiempo

#Modelos alfa generalizados
alpha_m = Constant(0.2)
alpha_f = Constant(0.4)
gamma   = Constant(0.5+alpha_f-alpha_m)
beta    = Constant((gamma+0.5)**2/4.)

#Parámetros de tiempo
T       = 4.0
Nsteps  = 100
dt = Constant(T/Nsteps)


#Malla
mesh = BoxMesh(Point(0., 0., 0.), Point(1., 0.1, 0.04), 60, 10, 5)

#Subdominios
def left(x, on_boundary):
    return near(x[0], 0.) and on_boundary

def right(x, on_boundary):
    return near(x[0], 1.) and on_boundary


# Propiedades del material
rho = Constant(1.0)

#Coeficientes de amortiguamiento de Rayleigh
eta_m = Constant(0.01)
eta_k = Constant(0.01)


#Carga aplicada
p0 = 1.
cutoff_Tc = T/5
#NOTAR COMO SE DEFINE LA FUNCIÓN!
p = Expression(("0", "t <= tc ? p0*t/tc : 0", "0"), t=0, tc=cutoff_Tc, p0=p0, degree=0)

#Se define igual que siempre
V = VectorFunctionSpace(mesh, 'Lagrange', degree=1)

#Funciones Test y trial
du = TrialFunction(V)
w = TestFunction(V)
# Desplazamiento actual Desconocido!
u = Function(V, name="Desplazamiento")
#Campos del paso anterior (Desplazamiento, velocidad, aceleración)
u0 = Function(V)
v0 = Function(V)
a0 = Function(V)


boundary_subdomains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundary_subdomains.set_all(0)
force_boundary = AutoSubDomain(right)
force_boundary.mark(boundary_subdomains, 3)

#Define la medida de la integral de borde (contorno)
dss = ds(subdomain_data=boundary_subdomains)

#Condición de borde del costado izquierdo
zero = Constant((0.0, 0.0, 0.0))
bc = DirichletBC(V, zero, left)

#Tensor de tensiones
def sigma(r):
    return 2.0*mu*sym(grad(r)) + lmbda*tr(sym(grad(r)))*Identity(len(r))

#Masa
def m(u, w):
    return rho*inner(u, w)*dx

#Rigidez
def k(u, w):
    return inner(sigma(u), sym(grad(w)))*dx

#Amortiguamiento de Rayleigh
def c(u, w):
    return eta_m*m(u, w) + eta_k*k(u, w)

#Trabajo de fuerzas externas
def Wext(w):
    return dot(w, p)*dss(3)

#---------------------------
#Ecuaciones de actualización
#CUIDADO!: pueden recibir tanto UFL como float

#Actualiza la aceleración
# a = 1/(2*beta)*((u - u0 - v0*dt)/(0.5*dt*dt) - (1-2*beta)*a0)
def actualiza_a(u, u0, v0, a0, ufl=True):
    if ufl:
        dt_ = dt
        beta_ = beta
    else:
        dt_ = float(dt)
        beta_ = float(beta)
    return (u-u0-dt_*v0)/beta_/dt_**2 - (1-2*beta_)/2/beta_*a0

#Acutaliza la velocidad
# v = dt * ((1-gamma)*a0 + gamma*a) + v0
def actualiza_v(a, u0, v0, a0, ufl=True):
    if ufl:
        dt_ = dt
        gamma_ = gamma
    else:
        dt_ = float(dt)
        gamma_ = float(gamma)
    return v0 + dt_*((1-gamma_)*a0 + gamma_*a)

def actualiza_campos(u, u0, v0, a0):
    """Actualiza los campos al final de cada paso.""" 

    # Get vectors (references)
    u_vec, u0_vec  = u.vector(), u0.vector()
    v0_vec, a0_vec = v0.vector(), a0.vector() 

    # use update functions using vector arguments
    a_vec = actualiza_a(u_vec, u0_vec, v0_vec, a0_vec, ufl=False)
    v_vec = actualiza_v(a_vec, u0_vec, v0_vec, a0_vec, ufl=False)

    # Update (u0 <- u)
    v0.vector()[:], a0.vector()[:] = v_vec, a_vec
    u0.vector()[:] = u.vector()



#Tiempos intermedio entre tn y tn+1 utilizando modelo alfa generalizado
def avg(x_old, x_new, alpha):
    return alpha*x_old + (1-alpha)*x_new

#Formulación variacional
a_new = actualiza_a(du, u0, v0, a0, ufl=True)
v_new = actualiza_v(a_new, u0, v0, a0, ufl=True)
res = m(avg(a0, a_new, alpha_m), w) + c(avg(v0, v_new, alpha_f), w) \
       + k(avg(u0, du, alpha_f), w) - Wext(w)
a_form = lhs(res)
L_form = rhs(res)


#Solvers (veremos luego)
K, res = assemble_system(a_form, L_form, bc)
solver = LUSolver(K, "mumps")
solver.parameters["symmetric"] = True


# Pasos de tiempo
time = np.linspace(0, T, Nsteps+1)
u_tip = np.zeros((Nsteps+1,))
energies = np.zeros((Nsteps+1, 4))
E_damp = 0
E_ext = 0
#Debemos definir un espacio de tensores extra:
Vsig = TensorFunctionSpace(mesh, "DG", 0)
sig = Function(Vsig, name="sigma")
xdmf_file = XDMFFile("resultados-elastodinamica.xdmf")
xdmf_file.parameters["flush_output"] = True
xdmf_file.parameters["functions_share_mesh"] = True
xdmf_file.parameters["rewrite_function_mesh"] = False


#Proyección para calcular la tensión
def local_project(v, V, u=None):
    """Element-wise projection using LocalSolver"""
    dv = TrialFunction(V)
    v_ = TestFunction(V)
    a_proj = inner(dv, v_)*dx
    b_proj = inner(v, v_)*dx
    solver = LocalSolver(a_proj, b_proj)
    solver.factorize()
    if u is None:
        u = Function(V)
        solver.solve_local_rhs(u)
        return u
    else:
        solver.solve_local_rhs(u)
        return

#iteración de pasos de tiempo
for (i, dt) in enumerate(np.diff(time)):

    t = time[i+1]
    print("Time: ", t)
    # Forces are evaluated at t_{n+1-alpha_f}=t_{n+1}-alpha_f*dt
    p.t = t-float(alpha_f*dt)

    # Solve for new displacement
    res = assemble(L_form)
    bc.apply(res)
    solver.solve(K, u.vector(), res)


    # Actualiza los campos, pero con valores numéricos
    actualiza_campos(u, u0, v0, a0)

    # Save solution to XDMF format
    xdmf_file.write(u, t)

    # Compute stresses and save to file
    local_project(sigma(u), Vsig, sig)
    xdmf_file.write(sig, t)

    p.t = t
    #Calculo de energias
    if MPI.comm_world.size == 1:
        u_tip[i+1] = u(1., 0.05, 0.)[1]
        
    E_elas = assemble(0.5*k(u0, u0))
    E_kin = assemble(0.5*m(v0, v0))
    E_damp += dt*assemble(c(v0, v0))
    E_tot = E_elas+E_kin+E_damp
    energies[i+1, :] = np.array([E_elas, E_kin, E_damp, E_tot])



np.savetxt('viscosidadM'+str(0.01)+'viscosidadK'+str(0.01)+'_dt'+str(dt)+'.out', (time, u_tip)) 

plt.figure(1)
plt.plot(time, u_tip,'r' ,label=r'$\alpha$-generalizado')
tiempo1_dt4,u_tip1_dt4 = np.loadtxt('viscosidad100.0_dt0.02.out')
plt.plot(tiempo1_dt4,u_tip1_dt4,'b',label= 'Tutorial anterior')
plt.xlabel("Tiempo")
plt.ylabel("Desplazamiento")
plt.legend(loc='upper right')
plt.ylim(-0.5, 0.5)



plt.figure(2)
plt.plot(time, energies)
plt.legend(("Elástica", "Cinética", "Amortiguamiento", "Total"))
plt.xlabel("Tiempo")
plt.ylabel("Energias")
plt.ylim(0, 0.0011)
plt.show()

