#!/usr/bin/env python
# coding: utf-8

#ejemplo23.py tomado de:
#https://comet-fenics.readthedocs.io/en/latest/demo/thermoelasticity/thermoelasticity_transient.html

from dolfin import *
#from mshr import *
import numpy as np
import matplotlib.pyplot as plt

L, H = 5, 0.3

mesh = RectangleMesh(Point(0., 0.), Point(L, H), 100, 10, "crossed")

T0 = Constant(293.)
DeltaT = Constant(10.)

E = Constant(50e3)
nu = Constant(0.2)
alpha = Constant(1e-5)
lmbda = Constant(E*nu/((1+nu)*(1-2*nu)))
mu = Constant(E/2/(1+nu))
rho = Constant(2700.)     # Densidad

kappa = Constant(alpha*(2*mu + 3*lmbda))
cV = Constant(910e-6)*rho  # Calor específico por unidad de volumen a deformación constante.
k = Constant(237e-6)  # Conductividad térmica

rho_g = 10e-3
b = Constant((0, -rho_g))

Vue = VectorElement('CG', mesh.ufl_cell(), 2) # elementos para el desplazamiento
Vte = FiniteElement('CG', mesh.ufl_cell(), 1) # elmentos para la temperatura
V = FunctionSpace(mesh, Vue*Vte)

def laterales(x, on_boundary):
    return (near(x[0], 0) or near(x[0], L)) and on_boundary
def abajo(x, on_boundary):
    return near(x[1], 0) and on_boundary
def arriba(x, on_boundary):
    return near(x[1], H) and on_boundary

bc1 = DirichletBC(V.sub(0), Constant((0., 0.)), laterales)
bc2 = DirichletBC(V.sub(1), Constant(0.), arriba)
bc3 = DirichletBC(V.sub(1), DeltaT, abajo)
bc4 = DirichletBC(V.sub(1), Constant(0.), laterales)
bcs = [bc1, bc2, bc3, bc4]

U_ = TestFunction(V)
(w_, v_) = split(U_) 
dU = TrialFunction(V)
(dw, dv) = split(dU) 
Uold = Function(V)
(uold, Thetaold) = split(Uold)

def eps(w):
    return sym(grad(w))


def sigma(w, Theta):
    return (lmbda*tr(eps(w)) - kappa*Theta)*Identity(2) + 2*mu*eps(w)

dt = Constant(0.)
formaVar_mecanica = inner(sigma(dw, dv), eps(w_))*dx-inner(b, w_)*dx
formVar_termica = (cV*(dv-Thetaold)/dt*v_ +
              kappa*T0*tr(eps(dw-uold))/dt*v_ +
              dot(k*grad(dv), grad(v_)))*dx

form = formaVar_mecanica + formVar_termica


Nincr = 100
t = np.logspace(1, 4, Nincr+1)
Nx = 100
x = np.linspace(0, L, Nx)
T_res = np.zeros((Nx, Nincr+1))
U = Function(V)


for (i, dti) in enumerate(np.diff(t)):
    print("Incremento: " + str(i+1))
    dt.assign(dti)
    solve(lhs(form) == rhs(form), U, bcs)
    Uold.assign(U)
    
    print('Temperatura en el centro: ',U(L/2, H/2)[2])
    
    T_res[:, i+1] = [U(xi, H/2)[2] for xi in x]


plt.figure()
plt.plot(x, T_res[:, 1::Nincr//10])
plt.xlabel("Coordenada $x$ en $y=H/2$")
plt.ylabel("Variación de la temperatura $\Theta$")
plt.legend(["$t={:.0f}$".format(ti) for ti in t[1::Nincr//10]], ncol=2)
plt.show()

u, Theta = split(U)

## Exportando los valores
#file_results = XDMFFile("datos_output.xdmf")
#file_results.parameters["flush_output"] = True
#file_results.parameters["functions_share_mesh"] = True

#file_results.write(u)
#file_results.write(sigma(u,Theta))
#file_results.write(Theta)

plt.figure()
p = plot(sigma(u, Theta)[1, 1], title="$\sigma_{yy}$")
#plt.xlim((0, 3*R))
#plt.ylim((0, 3*R))
plt.colorbar(p)
#plt.show()

plt.figure()
p = plot(1e3*u,title="Desplazamiento [mm]", mode="displacement")
plt.colorbar(p)

plt.figure()
p = plot(Theta, title="Variación de la temperatura")
#plt.xlim((0, 3*R))
#plt.ylim((0, 3*R))
plt.colorbar(p)



plt.show()
