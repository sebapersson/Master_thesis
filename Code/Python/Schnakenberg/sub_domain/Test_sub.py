#!/usr/bin/env python

from mshr import *
import numpy as np
import matplotlib.pyplot as plt
import meshio

# Convert the mesh to xdmf format 
msh = meshio.read("sub_dom.msh")
# Write all triangles
meshio.write("mesh.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]}))
# Write the triangle subdomains 
meshio.write("subdomains.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]},
                                    cell_data={"triangle": {"name_to_read": msh.cell_data["triangle"]["gmsh:physical"]}}))
# Write the lines 
meshio.write("lines.xdmf", meshio.Mesh(points=msh.points, cells={"line": msh.cells["line"]},
                                    cell_data={"line": {"name_to_read": msh.cell_data["line"]["gmsh:physical"]}}))


# Place dolfin here to avoid over-writing 
from dolfin import *

# Reading the mesh 
mesh = Mesh()
with XDMFFile("mesh.xdmf") as infile:
    infile.read(mesh)

# Reading the triangle subdomains and storing them 
sub_domains_pre = MeshValueCollection("size_t", mesh, 2)
with XDMFFile("subdomains.xdmf") as infile:
    infile.read(sub_domains_pre, "name_to_read")
sub_domains = cpp.mesh.MeshFunctionSizet(mesh, sub_domains_pre)

# Reading the line boundaries and storing them 
line_domain_pre = MeshValueCollection("size_t", mesh, 1)
with XDMFFile("lines.xdmf") as infile:
    infile.read(line_domain_pre, "name_to_read")
line_domain = cpp.mesh.MeshFunctionSizet(mesh, line_domain_pre)

# ----------------------------------------------------------------------------------
# Solving the problem in question
# ----------------------------------------------------------------------------------
# Class to hold the parameters for the Schankenberg model 
class param_schankenberg:
    def __init__(self,a=0.2,b=2.0,gamma=100.0,d=50.0):
        self.a = a
        self.b = b
        self.gamma = gamma
        self.d = d

# Class for constructing the initial conditions for each state (python syntax)
class ic(UserExpression):
    def __init__(self, *args, **kwargs):
        self.u0_1 = kwargs.pop('u0_1')
        self.u0_2 = kwargs.pop('u0_2')
        self.sigma = kwargs.pop('sigma')
        super(ic, self).__init__(*args, **kwargs)
        
    def eval(self,values,x):
        if (x[0]**2 + x[1]**2) < 0.5 + 1E-5:
            values[0] = 0.0
            values[1] = 0.0
        else:
            values[0] = self.u0_1 + np.random.normal(scale = self.sigma)
            values[1] = self.u0_2 + np.random.normal(scale = self.sigma)
        
    def value_shape(self):
        return(2,)


param = param_schankenberg(gamma = 10, d=100)
t_end = 2
n_time_step = 150

set_log_level(40)

# Initial values
u0_1 = param.a + param.b
u0_2 = (u0_1 - param.a)/u0_1**2

# Step size in t 
dt_inv = 1 / (t_end / n_time_step)
dt = t_end / n_time_step

# Standard Lagrange elements
P1 = FiniteElement('P', mesh.ufl_cell(), 1)
element = MixedElement([P1, P1])

# Function space to approximate solution in 
V = FunctionSpace(mesh, element)
u = Function(V)
u_n = Function(V)

# The initial conditions
u0_exp = ic(u0_1 = u0_1, u0_2 = u0_2, sigma = 0.05, element = V.ufl_element())

# Test and trial functions 
v_1, v_2 = TestFunctions(V)
u_1, u_2 = split(u)
u_n.interpolate(u0_exp)
u_n1, u_n2 = split(u_n)

# Define the different measures 
dx = Measure("dx", domain=mesh, subdomain_data=sub_domains)
ds = Measure("ds", domain=mesh, subdomain_data=line_domain)
dS = Measure("dS", domain=mesh, subdomain_data=line_domain)

# Check that the measures are correct
area_rec_not_circle = Constant(1.0) * dx(1)
area_circle = Constant(1.0) * dx(2)
c_circle = Constant(1.0) * dS(10)
area_rec_not_circle_ass = assemble(area_rec_not_circle)
area_circle_ass = assemble(area_circle)
c_circle_ass = assemble(c_circle)

print("Circle area = {:.3f}".format(area_circle_ass))
print("Rectangle without circle area = {:.3f}".format(area_rec_not_circle_ass))
print("Circle c = {:.3f}".format(c_circle_ass))

# The weak form, backward Euler in time  
dt_inv = Constant(dt_inv)
a = Constant(param.a)
b = Constant(param.b)
d = Constant(param.d)
gamma = Constant(param.gamma)

g = Constant(1E-100)

F = dt_inv*(u_1 - u_n1)*v_1*dx(1) - gamma*(a - u_1 + (u_1**2)*u_2)*v_1*dx(1) \
    + inner(grad(u_1), grad(v_1))*dx(1) + dt_inv*(u_2 - u_n2)*v_2*dx(1) \
    - gamma*(b - (u_1**2)*u_2)*v_2*dx(1) + d*inner(grad(u_2), grad(v_2))*dx(1)

#F += dt_inv*(u_1 - u_n1)*v_1*dx(2) - gamma*(a - u_1 + (u_1**2)*u_2)*v_1*dx(2) \
#    + inner(grad(u_1), grad(v_1))*dx(2) + dt_inv*(u_2 - u_n2)*v_2*dx(2) \
#    - gamma*(b - (u_1**2)*u_2)*v_2*dx(2) + d*inner(grad(u_2), grad(v_2))*dx(2)

#F += g("+") * v_1("+") * dS(10)
#F += g("+") * v_2("+") * dS(10)

# Files for saving the end result
folder_save = "test_sub_save"
file1 = folder_save + "/u_1.pvd"
file2 = folder_save + "/u_2.pvd"
vtkfile_u_1 = File(file1)
vtkfile_u_2 = File(file2)

# Solving the problem in time
t = 0
print("Starting to solve the PDE")
dU = Function(V)
for n in range(n_time_step):
    if (n + 1) % 10 == 0:
        print("Time step {} of {}".format(n+1, n_time_step))
        
    
    t += dt
    # We iterate 10 times. In future we might add a termination criteria
    for m in range(0, 10):
        A = assemble(derivative(F, u), keep_diagonal = True) # Define the derivative
        R = assemble(F) # Define the right hand side
        A.ident_zeros() # Fix "zero"-rows to get consistent system
        solve(A, dU.vector(), R) # Solver matrix system
        u.assign(u - dU) # Update solution one step    
    
    # Solve variational problem for time step
    # solve(F == 0, u)
    # Save solution to file 
    _u_1, _u_2, = u.split()
    vtkfile_u_1 << (_u_1, t)
    vtkfile_u_2 << (_u_2, t)
    u_n.assign(u)
