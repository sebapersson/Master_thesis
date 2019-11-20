#!/usr/bin/env python

from fenics import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt
import meshio

'''
This file aims to solve the Schankenberg reaction diffusion system using the finite 
element method (in space) and backward Euler in time. 
'''

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
        values[0] = self.u0_1 + np.random.normal(scale = self.sigma)
        values[1] = self.u0_2 + np.random.normal(scale = self.sigma)
        
    def value_shape(self):
        return(2,)

# ---------------------------------------------------------------------------------
# Start of functions 
# ---------------------------------------------------------------------------------

# Function that will solve the Schankenberg reaction diffusion system
# using the finite element method with Lagrange elements for the space
# parameters and a backward Euler for stepping in time
# Args:
#    param, the parameters for the Schnakenberg model (a, b, d, gamma)
#    t_end, the end time
#    n_time_step, the number of time steps
#    mesh, the mesh which the problem will be solved over
#    folder_save, the folder where the solution will be stored 
# Returns:
#    u, the values at the last time step (can be used to visualize the final time-step)
def solve_schnakenberg(param, t_end, n_time_step, mesh, folder_save):
    
    # Only display potential errors in the code
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
    
    # The weak form, backward Euler in time  
    dt_inv = Constant(dt_inv)
    a = Constant(param.a)
    b = Constant(param.b)
    d = Constant(param.d)
    gamma = Constant(param.gamma)
    F = dt_inv*(u_1 - u_n1)*v_1*dx - gamma*(a - u_1 + (u_1**2)*u_2)*v_1*dx \
        + inner(grad(u_1), grad(v_1))*dx + dt_inv*(u_2 - u_n2)*v_2*dx \
        - gamma*(b - (u_1**2)*u_2)*v_2*dx + d*inner(grad(u_2), grad(v_2))*dx
    
    # Files for saving the end result
    file1 = folder_save + "/u_1.pvd"
    file2 = folder_save + "/u_2.pvd"
    vtkfile_u_1 = File(file1)
    vtkfile_u_2 = File(file2)
    
    # Solving the problem in time
    t = 0
    print("Starting to solve the PDE")
    for n in range(n_time_step):
        if (n + 1) % 10 == 0:
            print("Time step {} of {}".format(n+1, n_time_step))
        
        t += dt
        # Solve variational problem for time step
        solve(F == 0, u)
        # Save solution to file 
        _u_1, _u_2, = u.split()
        vtkfile_u_1 << (_u_1, t)
        vtkfile_u_2 << (_u_2, t)
        u_n.assign(u)
    
    return u

# ---------------------------------------------------------------------------------
# End of functions, start of examples  
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# Rectangular mesh example
# ---------------------------------------------------------------------------------
param = param_schankenberg(gamma = 10, d=100)
t_end = 3
n_time_step = 150

# Rectangular mesh 
mesh = RectangleMesh(Point(0,0), Point(3, 3), 40, 40, "right/left")
path_save = "System_holes/reaction_system"
u = solve_schnakenberg(param, t_end, n_time_step, mesh, path_save)

# ---------------------------------------------------------------------------------
# Circular mesh example 
# ---------------------------------------------------------------------------------
param = param_schankenberg(d=100, gamma = 10)
t_end = 2
n_time_step = 150

# Create a circular mesh 
domain = Circle(Point(0, 0), 3)
mesh = generate_mesh(domain, 50)
path_save = "System_holes/reaction_system_circle"
u = solve_schnakenberg(param, t_end, n_time_step, mesh, path_save)

# ---------------------------------------------------------------------------------
# Rectangle with holes example 
# ---------------------------------------------------------------------------------
param = param_schankenberg(d=100, gamma=10)
t_end = 3
n_time_step = 150

box = Rectangle(Point(0, 0), Point(3, 3))
circle1 = Circle(Point(0.2, 0.2), 0.1)
circle2 = Circle(Point(0.3, 0.4), 0.1)
circle3 = Circle(Point(1, 1), 0.3)
domain = box - (circle1 + circle2 + circle3)
mesh = generate_mesh(domain, 50)
path_save = "System_holes/reaction_system_hole"
u = solve_schnakenberg(param, t_end, n_time_step, mesh, path_save)