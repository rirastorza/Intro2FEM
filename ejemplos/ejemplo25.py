#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

RFA_script: script that solve the bioheat equation

"""
from __future__ import print_function

import numpy as np
import time as tm
from scipy import constants as S
from fenics import *
import sys
import os
import matplotlib.pyplot as plt

import ufl as uflcond

pi = S.pi

# - Electric parameters
class EM_parameters:
    cond_rate = 0.
    V0 = 0. # ground pad voltage
    Vprobe = 100. # probe voltage
    conductivities = [0.28, 0.2, 1e-5 , 0.667]#Electric baseline conductivities of: Atrial wall, Connective tissue,Plastic, and Blood.
    Vmax = 66.33# Maximum voltage of the generator
    Vmin = 12.85# Minimum voltage of the generator

# - Thermal parameters
class thermal_parameters:
    kappas = [0.56,0.39,0.026,0.54]#  baseline thermal conductivities of: Atrial wall, Connective tissue,Plastic, and Blood.   
    rhoxcalorespe = [1081.*3686.,1027.*2372.,70.*1045.,1000.*4148.]  #Density times specific heat of: Atrial wall, Connective tissue,Plastic, and Blood.   
    dk = 0. # rate of change of thermal conductivity
    omega = Constant(0.004) # blood perfusion
    rho = Constant(1020.) # density blood
    c = Constant(3640.) # specific heat blood
    T0 = Constant(310.15) # baseline temperature
    T_initial = Constant(310.15) # initial flat temperature profile
    perf_model = 'stop'#
    cda_update = True # compute cell death
    p_stop = 1. # value of viability at which to stop perfusion
    ht = 610.0 # convective constant of the interface tissue-blood
    hp = 3446.0 #convective parameter of the plastic
    #Phase change parameters
    rho_c_t = 1060.*3411. # rho*c in tissue (phase)
    rho_c_v = 4.4e5 # rho*c in vapourised tissue (phase)
    Lh = 0 # latent heat of vapourisation
    Cliq = 0. # water content tissue (%)
    Tu = 374.-310.15 # upper transition temp
    Tl = 372.-310.15 # lower transition temp

# - Parameters of cell death Arrhenius function
class cell_death_parameters:
    A = 7.39e39
    R = S.R
    deltaE = 2.577e5

#Class to define the materials
class MATERIAL(UserExpression):
        def __init__(self, subdomains,s_1,s_2,s_3, s_4,**kwargs):
            super().__init__(**kwargs)
            self.subdomains = subdomains
            self.s_1 = s_1
            self.s_2 = s_2        
            self.s_3 = s_3
            self.s_4 = s_4        
        def eval_cell(self, values, x, cell):
            if self.subdomains[cell.index] == 1:
                values[0] = self.s_1#Tissue
            elif self.subdomains[cell.index] == 2:
                values[0] = self.s_2#Connective tissue
            elif self.subdomains[cell.index] == 3:
                values[0] = self.s_3#Plastic
            else:
                values[0] = self.s_4#Blood
                
def project_axisym(func,space):
    r = Expression('x[0]',degree=2)  
    w = TrialFunction(space)
    v = TestFunction(space)
    a = inner(w,v)*r*dx
    L = inner(func, v)*r*dx
    pfunc = Function(space)
    solve(a == L, pfunc)
    return pfunc 


#Function to compute de electric problem
def electric_problem(problemname, mesh, interior, boundaries, emp, Theta, thp):
    
    print("--+--+-- compute RFA SAR --+--+--")
    
    
    
    W = FunctionSpace(mesh, 'CG', 2)
    W_dg = FunctionSpace(mesh, 'DG', 2)
    W_dg0 = FunctionSpace(mesh, 'DG', 0) 

    T_p = project_axisym(Theta,W_dg0)# interpolate T onto piecewise constant to stabilise conductivity
    r = Expression('x[0]',degree=2)#Define the radius for axisymetric problem
    dss = ds(subdomain_data=boundaries)    # set measure
        
    # symmetry and insulating are natural BCs
    bcs = []
    bcs.append(DirichletBC(W, emp.Vprobe, boundaries, 10))#Active electrode
    bcs.append(DirichletBC(W, emp.Vprobe, boundaries, 40))#Active electrode at 40 Deg
    bcs.append(DirichletBC(W, emp.V0, boundaries, 20))#Passive electrode

    # define variational problem
    U = TrialFunction(W)
    V = TestFunction(W)
   
    # define electrical conductivity depending on temperature
    dependenceTsigma = (1.0 + emp.cond_rate*(T_p))
    
    #Diferent baseline conductivities
    conductivities = emp.conductivities
    s_1 = conductivities[0] #Atrial wall
    s_2 = conductivities[1] #Connective tissue
    s_3 = conductivities[2] #Plastic
    s_4 = conductivities[3] #Blood
    
    SIGMA = MATERIAL(interior,s_1,s_2,s_3, s_4,degree=0)
    
    a = dependenceTsigma*SIGMA*dot(grad(U), grad(V))*r*dx
    L = Constant(0.0)*V*r*dx

    U = Function(W)

    solve(a == L, U, bcs,solver_parameters={'linear_solver':'gmres','preconditioner':'ilu'})


    v = TestFunction(W_dg)
    qRF = TrialFunction(W_dg)
    a = qRF*v*r*dx
    L = v*dependenceTsigma*SIGMA*inner(nabla_grad(U), nabla_grad(U))*r*dx
    qRF = Function(W_dg)
    solve(a == L, qRF,solver_parameters={'linear_solver':'gmres','preconditioner':'ilu'})

    power = assemble(qRF*r*dx)*2.*np.pi
    resistance = (emp.Vprobe)**2/power
    
    print('Resistance: ',resistance)
    print('Voltage: ',emp.Vprobe)
    print('Power: ',power)

    return qRF, resistance, power, U

def rfa_bioheat_problem(mesh, interior, boundaries, problemname, dt,tmax, thp, emp):
    
    eps = np.finfo(float).eps 
    
    print("--+--+-- Start solving the bioheat equation --+--+--")
        
    W = FunctionSpace(mesh, 'CG', 2)
    W_dg = FunctionSpace(mesh, 'DG', 2) 
    W_dg_0 = FunctionSpace(mesh, 'DG', 0) 
    
    r = Expression('x[0]',degree=2)#Define the radius for axisymetric problem
    
    dss = ds(subdomain_data=boundaries)    # set measure
    
    # define quantities that need updating
    dte = Expression('dt', degree=1, dt=0.)
    cur_time = Expression('t', degree=1,t=0.)

    # initial uniform temperature
    Theta_prev = project_axisym(thp.T_initial-Constant(310.15),W)
    # initial values (if needed)
    resistance = 0.
        
    Q, resistance, power, Vfield = electric_problem(problemname, mesh, interior, boundaries, emp, Theta_prev, thp) 
    qext = project_axisym(Q,W)

    #Apply boundary conditions according to mesh function
    bcs = []
    bcs.append(DirichletBC(W, 0., boundaries, 20))
    bcs.append(DirichletBC(W, 0., boundaries, 50))
    bcs.append(DirichletBC(W, Constant(313.15)-thp.T_initial, boundaries, 40))#Constant temperature (40 Deg)
    
    #Define variables for variational form
    Theta_ = TrialFunction(W)
    v = TestFunction(W)
    f = qext 
       
    #Definition of thermal conductivity
    dependenceTkappa = (1.0 + thp.dk*(Theta_prev))    
    kappas = thp.kappas
    s_1 = kappas[0] #Atrial wall
    s_2 = kappas[1] #Connective tissue
    s_3 = kappas[2] #Plastic
    s_4 = kappas[3] #Blood
    KAPPA = MATERIAL(interior,s_1,s_2,s_3, s_4,degree=0)
    
    File("outputs/thermalconductivity.pvd") << project_axisym(KAPPA,W_dg_0)
    
    #Definitio of specific heat
    rhoxcalorespe = thp.rhoxcalorespe
    s_1 = rhoxcalorespe[0] #Atrial wall
    s_2 = rhoxcalorespe[1] #Connective tissue
    s_3 = rhoxcalorespe[2] #Plastic
    s_4 = rhoxcalorespe[3] #Blood
    
    RHOxCALESP = MATERIAL(interior,s_1,s_2,s_3, s_4,degree=0)
    
    #Definition of perfusion
    D_prev=interpolate(Constant(0.),W)
    omega=thp.omega
    D_prev_const = project_axisym(D_prev,W_dg_0)# project D onto piecewise constant mesh to stop negative values
    if thp.perf_model=='stop':
        omega=conditional(gt(D_prev_const,thp.p_stop), thp.omega, 0.)
        print("check perfusion threshhold")
        
    
    rc1 = conditional(lt(Theta_prev,thp.Tl), thp.rho_c_t, 0.) 
    #Modificado segun:
    #https://fenicsproject.org/docs/ufl/1.5.0/ufl.html#module-ufl.conditional
    rc2 = conditional(uflcond.operators.And(ge(Theta_prev,thp.Tl),le(Theta_prev,thp.Tu)), (thp.rho_c_t+thp.rho_c_v)/2+thp.rho*thp.Lh*thp.Cliq*(1./(thp.Tu-thp.Tl)), 0.)
    rc3 = conditional(gt(Theta_prev,thp.Tu), thp.rho_c_v, 0.)
      
    
    # Heat transfer variational form
    a = KAPPA*dependenceTkappa*inner(nabla_grad(Theta_), nabla_grad(v))*r*dx+ v*omega*thp.rho*thp.c*Theta_*r*dx+thp.hp*Theta_*v*r*dss(30) + v*rc1/dte*Theta_*r*dx + v*rc2/dte*Theta_*r*dx + v*rc3/dte*Theta_*r*dx
    L = f*v*r*dx + v*rc1/dte*Theta_prev*r*dx + v*rc2/dte*Theta_prev*r*dx + v*rc3/dte*Theta_prev*r*dx
        
    Theta_ = Function(W)

    # assemble in advance of time iteration
    A = None
    b = None

    Q = Function(W_dg)
        
    store_resistance = []# save the resistance at output times
    store_power = []# save power
    store_sensortemp = []# guardar temperatura en sensor
    tiempo = []# save the time
    power = 0.
    
    # initialise cell death
    n = len(Theta_.vector().get_local())
    cda = np.zeros(n) # cell death array
    D = interpolate(Constant(0.),W) # dead field
    
    t = dt
    
    #vtkfile = File('data/solution.pvd')
    while t <= tmax+eps:
        print("--+--+--             update qRF              --+--+--")
        Q, resistance, power, Vfield = electric_problem(problemname, mesh, interior, boundaries, emp, Theta_prev, thp)
        # assemble each iteration to account for previous time step
        dte.dt = dt
        cur_time.t = t
        f.assign(Q)
        b = assemble(L, tensor=b)
        A = assemble(a, tensor=A)
        for bc in bcs:
            bc.apply(A, b)
        solve(A, Theta_.vector(), b)
        
        sensorTemperature = Theta_(0.001,-0.00116)# Punta del electrodo

        #vtkfile << (project_axisym(Theta_+Constant(310.),W), t)

        nodal_T = Theta_.vector().get_local()
        nodal_T_prev = Theta_prev.vector().get_local()
        T_error = np.abs(nodal_T-nodal_T_prev).max()
        
        store_resistance.append(resistance)
        store_power.append(power)
        store_sensortemp.append(sensorTemperature)
        tiempo.append(t)

        print("***************************** CELL DEATH  **********************************")
                
        if thp.cda_update:
            cda = cell_death_timestep(cda,n,t,dt,nodal_T,cell_death_parameters)
        
        D.vector()[:] = cda
        D_prev.assign(D) # update cell death for perfusion
        
        CDeath = project_axisym(D_prev,W)
        
        print("Difference Temperature step before: ",T_error,"   t: ", t+dt, "   dt: ", dt, "   pow: ", power, "   imp: ", resistance)
               
        
        print("***************************** PI Controller  **********************************")
        # -------------
        # PI controller
        Theta_target = 353.15-310.15 # Target temperature = 80°C 
        Vmax = emp.Vmax #
        Vmin = emp.Vmin
        Kp = 280.78 #Proportional constant of the controller
        Ki = 0.5 #Integration constant of the controller
        xik1 = 0.0 #initial value of the integrator
        ek = Theta_target-sensorTemperature
        xik = xik1+ek #integrator
        Pout = Kp*ek+Ki*xik # output of the controller
        
        if Pout < 0.0:
            Pout = 0.0
        else:
            print('-----------------')
        
        Vpi = Pout**0.5#Voltage of the controller
        
        if Vpi > Vmax:
            Vpi = Vmax
            xik = xik1 #not integrate if satured
        elif Vpi < Vmin:
            Vpi = Vmin
            xik = xik1 #not integrate if satured
        else:
            print('------------------')
            
        xik1 = xik
        
        emp.Vprobe = Vpi
                
        t += dt
        Theta_prev.assign(Theta_)
        
    return Theta_+310.15-273.15, project_axisym(Vfield,W),project_axisym(Theta_+310.15-273.15,W_dg_0), store_resistance,store_sensortemp,store_power,tiempo,CDeath



def cell_death_func_class(t,y,args): # ode class
    # scaled for T' = T-310.15
    # args = n, cdp.deltaE, cdp.R, cdp.A, T_prev_nodal
    dydt = np.array(y)
    dydt = args[3]*np.exp(-args[1]/((args[4]+310.15)*args[2]))
    
    return dydt


def cell_death_timestep(y0,n,t,dt,T_prev_nodal,cdp):
    import scipy.integrate as spint
    step_int = spint.ode(cell_death_func_class)
    step_int.set_integrator("vode",method="bdf",nsteps=1e5)
    step_int.set_f_params([n, cdp.deltaE, cdp.R, cdp.A, T_prev_nodal])
    step_int.set_initial_value(y0,0.)
    step_int.integrate(dt)
    print(step_int.successful())
    if not step_int.successful():
        error("cell death step solve failed")
    cda = step_int.y
    return cda

start_time = tm.strftime('%H:%M:%S')
problemname = "rfa-tissue"


##Convierto geometria creada con gmsh
#string = "gmsh -2 "+problemname+".geo "
#os.system(string)

#string = "dolfin-convert "+problemname+".msh "+problemname+".xml"
#os.system(string)



EM_parameters.V0 = 0.
EM_parameters.Vprobe = 20.0
EM_parameters.cond_rate = 0.015 #1.5%/°C



# set thermal parameters
thermal_parameters.dk = Constant(0.02)
thermal_parameters.omega = Constant(0.004)
thermal_parameters.rho = Constant(1020.)
thermal_parameters.c = Constant(3640.)
thermal_parameters.T0 = Constant(310.15)
thermal_parameters.T_initial = Constant(310.15)

# solver options
tmax = 5.0 # maximum time (s)
dt = .1 # time step (s)
thermal_parameters.perf_model = 'stop'

set_log_active(False) # switch off fenics messages

# Load geometry
mesh = Mesh(problemname+".xml")
boundaries = MeshFunction("size_t", mesh, problemname+"_facet_region.xml")
interior = MeshFunction("size_t", mesh, problemname+"_physical_region.xml")


T,V,Tdeg, impedancia1, temperaturaSensor1,potenciaSensor1,tiempo1,CDeath = rfa_bioheat_problem(mesh, interior, boundaries, problemname, dt,tmax,thermal_parameters, EM_parameters)

File("outputs/temperature.pvd") << Tdeg
File("outputs/voltage.pvd") << V
File("outputs/CellDeath.pvd") << CDeath

print('start time: ', start_time)
print('end time:   ', tm.strftime('%H:%M:%S'))



from matplotlib import pyplot as plt

fig1 = plt.figure(1)
f1 = fig1.add_subplot(211)
f1.plot(np.asarray(tiempo1),np.asarray(impedancia1),'.r')
f1.set_title('Impedance')
f1 = fig1.add_subplot(212)
f1.plot(np.asarray(tiempo1),np.asarray(temperaturaSensor1)+37.,'.r')
f1.plot(tiempo1,80.*np.ones_like(np.asarray(impedancia1)),'k')
f1.set_title('Temperature of the sensor')
f1.set_xlabel('Time (s)')

tiempo = np.asarray(tiempo1)
Pot = np.asarray(potenciaSensor1)
Tem = np.asarray(temperaturaSensor1)
Imp = np.asarray(impedancia1)
np.savetxt('outputs/sensorData.out', (tiempo,Pot,Tem,Imp)) 

plt.show()

plt.figure()
p = plot(Tdeg, title="Final temperature (Celsius)")

plt.colorbar(p)

plt.show()



