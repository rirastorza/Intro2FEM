"""
Ejemplo21.py modificado del tutorial de FEniCS:

Seguimos con fluido incompresible. La formulación dura es:

  u' + u . nabla(u)) - div(sigma(u, p)) = f
                                 div(u) = 0
                                 
FEniCS tutorial demo program: Incompressible Navier-Stokes equations
for flow around a cylinder using the Incremental Pressure Correction
Scheme (IPCS).

"""

from __future__ import print_function
from fenics import *
from mshr import *
import numpy as np

T = 5.0            # tiempo final
num_steps = 5000   # número de pasos
dt = T / num_steps # tamaño de paso de tiempo
mu = 0.001         # viscosidad
rho = 1            # densidad

# Creo la malla
channel = Rectangle(Point(0, 0), Point(2.2, 0.41))
cylinder = Circle(Point(0.2, 0.2), 0.05)
domain = channel - cylinder
mesh = generate_mesh(domain, 64)

#plot (mesh)



# Notar que uno es vectorial y el otro escalar
V = VectorFunctionSpace(mesh, 'P', 2)
Q = FunctionSpace(mesh, 'P', 1)

# Condiciones de borde
inflow   = 'near(x[0], 0)'
outflow  = 'near(x[0], 2.2)'
walls    = 'near(x[1], 0) || near(x[1], 0.41)'
cylinder = 'on_boundary && x[0]>0.1 && x[0]<0.3 && x[1]>0.1 && x[1]<0.3'

# Perfil de entrada
inflow_profile = ('4.0*1.5*x[1]*(0.41 - x[1]) / pow(0.41, 2)', '0')

# Defino las condiciones de borde
bcu_inflow = DirichletBC(V, Expression(inflow_profile, degree=2), inflow)
bcu_walls = DirichletBC(V, Constant((0, 0)), walls)
bcu_cylinder = DirichletBC(V, Constant((0, 0)), cylinder)
bcp_outflow = DirichletBC(Q, Constant(0), outflow)
bcu = [bcu_inflow, bcu_walls, bcu_cylinder]
bcp = [bcp_outflow]

u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

#Defino las funciones para las soluciones en el paso previo y el actual
u_n = Function(V)
u_  = Function(V)
p_n = Function(Q)
p_  = Function(Q)


U  = 0.5*(u_n + u)
n  = FacetNormal(mesh)
f  = Constant((0, 0))
k  = Constant(dt)
mu = Constant(mu)
rho = Constant(rho)

#Tensor de velocidad de deformación (simétrico)
def epsilon(u):
    return sym(nabla_grad(u))

#Tensor de tensiones
def sigma(u, p):
    return 2*mu*epsilon(u) - p*Identity(len(u))

# formulación variacional para el paso 1
F1 = rho*dot((u - u_n) / k, v)*dx \
   + rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx \
   + inner(sigma(U, p_n), epsilon(v))*dx \
   + dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds \
   - dot(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# formulación variacional para el paso 2
a2 = dot(nabla_grad(p), nabla_grad(q))*dx
L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx

# formulación variacional para el paso 3
a3 = dot(u, v)*dx
L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx


A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Aplico condiciones de borde
[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]

## Creo XDMF para visualización
#xdmffile_u = XDMFFile('navier_stokes_cylinder/velocity.xdmf')
#xdmffile_p = XDMFFile('navier_stokes_cylinder/pressure.xdmf')

pdvfile_p = File("data/presion.pvd")
pdvfile_u = File("data/velocidad.pvd")

##---------------------------------------
## Interesante! Para usar en otro sistema!
#timeseries_u = TimeSeries('data/velocity_series')
#timeseries_p = TimeSeries('data/pressure_series')
## Salvo también la malla
#File('data/cylinder.xml.gz') << mesh
##-----------------------------------------------------

# La progress bar!
progress = Progress('Time-stepping', num_steps)

import matplotlib.pyplot as plt


plt.figure(1)
plot(mesh)
plt.show()

t = 0


for n in range(num_steps):

    # Update current time
    t += dt

    # Paso 1: Tentative velocity step
    b1 = assemble(L1)
    [bc.apply(b1) for bc in bcu]
    solve(A1, u_.vector(), b1, 'bicgstab', 'hypre_amg')

    # Paso 2: Pressure correction step
    b2 = assemble(L2)
    [bc.apply(b2) for bc in bcp]
    solve(A2, p_.vector(), b2, 'bicgstab', 'hypre_amg')

    # Paso 3: Velocity correction step
    b3 = assemble(L3)
    solve(A3, u_.vector(), b3, 'cg', 'sor')

    ## Plot solution
    #plt.figure(2)
    #plot(u_, title='Velocity')
    #plt.figure(3)
    #plot(p_, title='Pressure')
    #plt.show()

    ## Save solution to file (XDMF/HDF5)
    #xdmffile_u.write(u_, t)
    #xdmffile_p.write(p_, t)
    
    pdvfile_p << p_
    pdvfile_u << u_
    
    ## Save nodal values to file
    #timeseries_u.store(u_.vector(), t)
    #timeseries_p.store(p_.vector(), t)

    # Update previous solution
    u_n.assign(u_)
    p_n.assign(p_)

    # Update progress bar
    set_log_level(LogLevel.PROGRESS)
    progress += 1
    set_log_level(LogLevel.ERROR)
    print('u max:', u_.vector().get_local().max())

# Hold plot
interactive()

list_krylov_solver_preconditioners()
