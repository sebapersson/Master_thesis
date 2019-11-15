#!/usr/bin/env python

from fenics import *
import numpy as np
import matplotlib.pyplot as plt

'''
This file aims to solve the Schankenberg reaction diffusion system using the finite 
element method (in space) and backward Euler in time. 
'''

# The model parameters 
gamma = 1000
d = 50
a = 0.2
b = 2

# The initial values
u0_1 = a + b
u0_2 = (u0_1 - a)/u0_1**2

#set_log_level(40)

# Solution parameters
T = 1
num_steps = 300
dt_inv = 1 / (T / num_steps)
dt = T / num_steps

# Define a triangular mesh on a rectangle 
mesh = RectangleMesh(Point(0,0), Point(3, 2), 30, 30, "right/left")

# Linear functions on a triangular mesh 
P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1])
V = FunctionSpace(mesh, element)

# Defining the test functions for the problem
v_1, v_2 = TestFunctions(V)

# The approximate functions
u = Function(V)
u_1, u_2 = split(u)

# The initial values
u_n = Function(V)
u0_exp = Expression(('u0_1 + 0.1*sin(x[0]) + 0.1*sin(x[1])', 'u0_2 + 0.01'),
                    u0_1=u0_1, u0_2=u0_2, element = element)
u_n = interpolate(u0_exp, V)
u_n1, u_n2 = split(u_n)

# Expressions used in the variational form
dt_inv = Constant(dt_inv)
a = Constant(a)
b = Constant(b)
d = Constant(d)
gamma = Constant(gamma)

# The finite element formulation using backward Euler method 
F = dt_inv*(u_1 - u_n1)*v_1*dx - gamma*(a - u_1 + (u_1**2)*u_2)*v_1*dx \
    + inner(grad(u_1), grad(v_1))*dx + dt_inv*(u_2 - u_n2)*v_2*dx \
    - gamma*(b - (u_1**2)*u_2)*v_2*dx + d*inner(grad(u_2), grad(v_2))*dx

# For saving the result
vtkfile_u_1 = File('reaction_system/u_1.pvd')
vtkfile_u_2 = File('reaction_system/u_2.pvd')

# Solving the problem in time
t = 0
for n in range(num_steps):
    # Update current time
    t += dt
    # Solve variational problem for time step
    solve(F == 0, u)
    # Save solution to file 
    _u_1, _u_2, = u.split()
    vtkfile_u_1 << (_u_1, t)
    vtkfile_u_2 << (_u_2, t)
    u_n.assign(u)
