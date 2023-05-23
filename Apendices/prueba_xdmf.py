from dolfin import *

nx = ny = 8
mesh = UnitSquareMesh(nx, ny)
V = FunctionSpace(mesh, 'P', 1)

alpha = 3; beta = 1.2
f = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t', degree=2, alpha=alpha, beta=beta, t=0)
f_out = XDMFFile("test.xdmf")

f_out.write_checkpoint(project(f,V), "f", 0.0,  XDMFFile.Encoding.HDF5, False)

for j in range(1,5):
    t = float(j/5)
    f.t = t
   
    f_out.write_checkpoint(project(f,V), "f", t,  XDMFFile.Encoding.HDF5, True)

f_out.close()

f1 = Function(V)
f_in =  XDMFFile("test.xdmf")

f_in.read_checkpoint(f1,"f",0) 
f_in.read_checkpoint(f1,"f",1) 
f_in.close()
