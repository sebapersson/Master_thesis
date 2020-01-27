#!/usr/bin/env python

import numpy as np
import pandas as pd
from dolfin import *
import os
import sys
from tqdm import tqdm
import matplotlib.pyplot as plt


# Class to hold the time options 
class t_opt_class:
    def __init__(self, t_end, n_time_step):
        self.t_end = t_end
        self.n_time_step = n_time_step

# Class to hold the parameters for the Gierer-Meinhardt model
class param_gierer:
    def __init__(self, a=0.5, b=2.0, gamma=100.0, d=50.0):
        self.a = a
        self.b = b
        self.gamma = gamma
        self.d = d

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
        self.x_max = kwargs.pop('x_max')
        super(ic, self).__init__(*args, **kwargs)
        
    def eval(self,values,x):
        if x[0] > (self.x_max / 2 - 0.1) and x[0] < (self.x_max / 2 + 0.1):
            values[0] = self.u0_1 + np.absolute(np.random.normal(scale = self.sigma))
            values[1] = self.u0_2
        else:
            values[0] = self.u0_1
            values[1] = self.u0_2
    
    def value_shape(self):
        return(2,)


# Function that will solve the Schankenberg system for a one-dimensional
# case, the purpose of this is to act as an illustration of how the
# reaction diffusion models work. # TODO: extend to Gierer model
# Args:
#    param , the model parameters
#    t_opt, a t-opt class object (contains time-data)
#    x_max, the maximum x-coordinate
#    dir_save, the path to the directory to where the result is saved
#    file_save, path to the file where result is saved
#    model, the model that is solved 
#    seed, the model seed (deafult 123)
def fem_illustration(param, t_opt, dir_save, file_save, model="Schankenberg", x_max=10, seed=123):
    
    np.random.seed(seed)
    t_end = t_opt.t_end
    n_time_step = t_opt.n_time_step
    
    # Time parameters 
    dt_inv = 1 / (t_end / n_time_step)
    dt = t_end / n_time_step
    
    # The initial steady state
    if model == "Schankenberg":
        u0_1 = param.a + param.b
        u0_2 = (u0_1 - param.a) / (u0_1**2)
    elif model == "Gierer":
        u0_1 = (param.a + 1) / param.b
        u0_2 = u0_1 ** 2
    else:
        print("Error: model {} isn't an option".format(model))
        sys.exit(1)
    
    # Function space to solve problem in 
    mesh = IntervalMesh(1500, 0, x_max)
    P1 = FiniteElement('P', mesh.ufl_cell(), 1)
    element = MixedElement([P1, P1])
    
    # Function space to approximate solution in
    V = FunctionSpace(mesh, element)
    
    # Expression for the initial conditions 
    u0_exp = ic(u0_1=u0_1, u0_2=u0_2, sigma=0.05, x_max=x_max, element=V.ufl_element())
    
    # Test functions, and initial values
    u_n = Function(V)
    v_1, v_2 = TestFunctions(V)
    u_n.interpolate(u0_exp)
    u_n1, u_n2 = split(u_n)
    
    # Formulate the weak-formulation
    u_1, u_2 = TrialFunctions(V)
    dt_inv = Constant(dt_inv)
    a = Constant(param.a)
    b = Constant(param.b)
    d = Constant(param.d)
    gamma = Constant(param.gamma)
    
    if model == "Schankenberg":
        F = dt_inv*(u_1 - u_n1)*v_1*dx - gamma*(a - u_n1 + (u_n1**2)*u_n2)*v_1*dx \
            + inner(grad(u_1), grad(v_1))*dx + dt_inv*(u_2 - u_n2)*v_2*dx \
            - gamma*(b - (u_n1**2)*u_n2)*v_2*dx + d*inner(grad(u_2), grad(v_2))*dx
    elif model == "Gierer":
        F = dt_inv*(u_1 - u_n1)*v_1*dx - gamma*(a - b * u_n1 + (u_n1**2) / u_n2)*v_1*dx \
            + inner(grad(u_1), grad(v_1))*dx + dt_inv*(u_2 - u_n2)*v_2*dx \
            - gamma*(u_n1 ** 2 - u_n2)*v_2*dx + d*inner(grad(u_2), grad(v_2))*dx
    
    t = 0
    U = Function(V)
    U.assign(u_n)
    
    # For saving the x-coordinates
    n = V.dim()  
    d = mesh.geometry().dim()
    dof_coordinates = V.tabulate_dof_coordinates().reshape(n, d)
    dof_x = dof_coordinates[:, 0]
    
    # Create file to where the result will be saved (remove if present)
    if not os.path.isdir(dir_save):
        os.mkdir(dir_save)
    if os.path.isfile(file_save):
        os.remove(file_save)
    
    # Solving a linear system, all non-linear into last-vector 
    FA = lhs(F)
    Fb = rhs(F)
    A = assemble(FA, keep_diagonal = True)
    A.ident_zeros()
    time_vec = np.linspace(0, t_end, n_time_step)
    
    # For writing to file (save initial state)
    _u_1, _u_2, = U.split()
    u1_vec = _u_1.vector().get_local()
    u2_vec = _u_2.vector().get_local()
    data_to_save = pd.DataFrame({"t": t, "u1": u1_vec, "u2": u2_vec, "x": dof_x})
    
    n = 0
    for n in tqdm(range(n_time_step)):
        
        # Update time-step
        t += dt
        
        # Solve linear variational problem for time step
        b = assemble(Fb)
        solve(A,  U.vector(), b)        
        
        _u_1, _u_2, = U.split()
        
        # Sanity check solution
        u1_vec = _u_1.vector().get_local()
        u2_vec = _u_2.vector().get_local()
        min_u1 = np.min(u1_vec)
        min_u2 = np.min(u2_vec)
        if min_u1 < 0 or min_u2 < 0:
            print("Error, negative concentration")
            sys.exit(1)
        
        # Update previous solution 
        u_n.assign(U)
        
        # For saving result
        if n % 20 == 0:
            df = pd.DataFrame({"t": t, "u1": u1_vec, "u2": u2_vec, "x": dof_x})
            data_to_save = data_to_save.append(df)
        n += 1
    
    # Write the result to file 
    n_rep = int(len(data_to_save.index) / 2)
    test = np.tile([1, 2], n_rep)
    data_to_save["id_mol"] = test
    data_to_save.to_csv(file_save)
    

# Schankenberg model, run illustration case 
t_opt = t_opt_class(7.5, 1500)
param = param_schankenberg(gamma=5, d=200)
dir_save = "../../../Intermediate/Illustration/"
file_save = dir_save + "Schankenberg_illustration.csv"
#fem_illustration(param, t_opt, dir_save, file_save)

# Gierer model, run illustration case 
t_opt = t_opt_class(1.5, 2000)
param = param_gierer(b = 2.0, a = 0.5, gamma = 20, d = 50)
dir_save = "../../../Intermediate/Illustration/"
file_save = dir_save + "Gierer_illustration.csv"
fem_illustration(param, t_opt, dir_save, file_save, x_max= 4, model="Gierer")

