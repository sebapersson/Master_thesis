#!/usr/bin/env python

import numpy as np
import pandas as pd
from dolfin import *
import meshio
import os
import sys
import time 
from tqdm import tqdm

# -----------------------------------------------------------------------------------
# Classes for solving the problem
# -----------------------------------------------------------------------------------

# Class to hold the information when sanity checking the solutions, tag is a user-added
# ending to know which experiment is sanity checked. 
class file_loc_sanity_check:
    def __init__(self, geom, n_holes, model, tag):
        if n_holes == "Zero_holes":
            hole = "h0_"
        elif n_holes == "Five_holes":
            hole = "h5_"
        elif n_holes == "Seven_holes":
            hole = "h7_"
        elif n_holes == "Twenty_holes":
            hole = "h20_"
        else:
            print("Error, not allowed number of holes")
            sys.exit(1)
        self.mesh_folder = "../../../Intermediate/" + geom + "_mesh/" + n_holes + "_mesh/"
        self.pwd_folder = "../../../Result/" + model + "/" + geom + "/pwd_files/sanity_check/" + tag + "/" + hole
        self.model = model

# Class object to hold which sub-domains will have changed parameters
# values, and what the parameter values are changed to. 
class diff_param_class:
    def __init__(self, dx_list, param):
        self.param = param
        self.dx_list = dx_list
    
    # Method to convert the parameters to string
    def convert_to_str(self):
        a = str(self.param.a).replace(".", "D")
        b = str(self.param.b).replace(".", "D")
        gamma = str(self.param.gamma).replace(".", "D")
        d = str(self.param.d).replace(".", "D")
        return a, b, gamma, d


# Class to hold the time options 
class t_opt_class:
    def __init__(self, t_end, n_time_step):
        self.t_end = t_end
        self.n_time_step = n_time_step


# Class for holding the initial value parameters, if geom is equal to rectangle
# a quadratic rectangle will be disturbed with middle x_mid, y_mid and side length r_circle*2
class init_val_param:
    def __init__(self, controlled=False, geom="Circle", x_mid=0, y_mid=0, r_circle=0.25):
        self.controlled = controlled
        self.x_mid = x_mid
        self.y_mid = y_mid
        self.r_circle = r_circle
        self.geom = geom
    
    # Method that converts and returns the arguments as strings
    def convert_str(self):
        r = str(self.r_circle).replace(".", "D")
        x = str(self.x_mid).replace(".", "D")
        y = str(self.y_mid).replace(".", "D")
        return r, x, y


# Class for constructing the initial conditions for each state (python syntax)
# Note that the pattern formation is sensitive to the noise sigma. 
class ic(UserExpression):
    def __init__(self, *args, **kwargs):
        self.u0_1 = kwargs.pop('u0_1')
        self.u0_2 = kwargs.pop('u0_2')
        self.sigma = kwargs.pop('sigma')
        self.init_par = kwargs.pop('init_param')
        super(ic, self).__init__(*args, **kwargs)
        
    def eval(self,values,x):
        r = self.init_par.r_circle
        x_mid = self.init_par.x_mid; y_mid = self.init_par.y_mid
        if self.init_par.geom == "Circle" and self.init_par.controlled == True:
            if (x[0] - x_mid)**2 + (x[1] - y_mid)**2 < r**2:
                values[0] = self.u0_1 + np.absolute(np.random.normal(scale = self.sigma))
                values[1] = self.u0_2 + np.absolute(np.random.normal(scale = self.sigma))
            else:
                values[0] = self.u0_1
                values[1] = self.u0_2
        # The case the disturbance is a rectangle 
        elif self.init_par.controlled == True and self.init_par.geom == "Rectangles":
            # Extreme values of the rectangle
            ext_val1 = -r
            ext_val2 = r
            if (x[0] > ext_val1 and x[0] < ext_val2) and (x[1] > ext_val1 and x[1] < ext_val2):
                values[0] = 1.1 * self.u0_1 
                values[1] = 1.1 * self.u0_2 
            else:
                values[0] = self.u0_1
                values[1] = self.u0_2
        # Case the entire region will be disturbed 
        else:
            values[0] = self.u0_1 + np.random.normal(scale = self.sigma)
            values[1] = self.u0_2 + np.random.normal(scale = self.sigma)
        
    def value_shape(self):
        return(2,)


# Class to hold the parameters for the Schankenberg model 
class param_schankenberg:
    def __init__(self, a=0.2, b=2.0, gamma=100.0, d=50.0):
        self.a = a
        self.b = b
        self.gamma = gamma
        self.d = d

# Class to hold the parameters for the Gierer-Meinhardt model
class param_gierer:
    def __init__(self, a=0.5, b=2.0, gamma=100.0, d=50.0):
        self.a = a
        self.b = b
        self.gamma = gamma
        self.d = d

# Class to hold the file-locations for a model, note that diff_para_original is a string that
# will contain information if any of the standard parameters have been changed
class file_locations_class:
    def __init__(self, n_holes, geometry, model, ic_par, controlled_inital=False, find_mesh_size=False,
                 lc="0", diff_para=None):
        # If the initial values aren't controlled and not finding mesh-size
        if n_holes == "Zero_holes":
            hole = "0"
        elif n_holes == "Five_holes":
            hole = "5"
        elif n_holes == "Seven_holes":
            hole = "7"
        elif n_holes == "Twenty_holes":
            hole = "20"
        else:
            print("Error: Invalid number of holes")
            print("n_holes = {}".format(n_holes))
            sys.exit(1)
        
        # Information about disturbing initial values
        r, x_mid, y_mid = ic_par.convert_str()
        
        if controlled_inital == False and find_mesh_size == False and diff_para == None:
            self.path_to_msh_file = "../../Gmsh/" + geometry + "/" + n_holes + ".msh"
            self.mesh_folder = "../../../Intermediate/" + geometry + "_mesh/" + n_holes + "_mesh/"
            self.pwd_folder = "../../../Result/" + model + "/" + geometry + "/pwd_files/" + "h" + hole + "_d_k"
            self.file_save_folder = "../../../Intermediate/" + model + "_files/" + geometry + "/h" + hole + "_d_k"
            self.model = model
        # Rename save-folders if initial values are controlled 
        elif controlled_inital == True and find_mesh_size == False and diff_para == None:
            self.path_to_msh_file = "../../Gmsh/" + geometry + "/" + n_holes + ".msh"
            self.mesh_folder = "../../../Intermediate/" + geometry + "_mesh/" + n_holes + "_mesh/"
            self.pwd_folder = ("../../../Result/" + model + "/" + geometry + "/pwd_files/" + "h" + hole + 
                              "_d" + "r" + r + "x" + x_mid + "y" + y_mid +  "_k")
            self.file_save_folder = ("../../../Intermediate/" + model + "_files/" + geometry + "/h" + hole + 
                              "_d" + "r" + r + "x" + x_mid + "y" + y_mid +  "_k")
            self.model = model
        # Rename save-folders if the aim is to find the mesh-size
        elif find_mesh_size == True:
            # Change . to D structure
            lc = lc.replace('.', 'D')
            # Create the files 
            self.path_to_msh_file = "../../Gmsh/" + geometry + "/" + n_holes + ".msh"
            self.mesh_folder = "../../../Intermediate/" + geometry + "_mesh/" + n_holes + "_mesh/"
            self.pwd_folder = "../../../Result/" + model + "/" + geometry + "/pwd_files/" + "h" + hole + "_lc" + lc 
            self.file_save_folder = "../../../Intermediate/" + model + "_files/" + geometry + "/h" + hole + "_lc" + lc 
            self.model = model
        # The case the parameters have been disturbed
        elif diff_para != None and controlled_inital == False:
            a, b, gamma, d = diff_para.convert_to_str()
            self.path_to_msh_file = "../../Gmsh/" + geometry + "/" + n_holes + ".msh"
            self.mesh_folder = "../../../Intermediate/" + geometry + "_mesh/" + n_holes + "_mesh/"
            self.pwd_folder = ("../../../Result/" + model + "/" + geometry + "/pwd_files/" + "h" + hole + "_d_k" +
                               "a" + a + "b" + b + "ga" + gamma + "di" + d)
            self.file_save_folder = ("../../../Intermediate/" + model + "_files/" + geometry + "/h" + hole + "_d_k" +
                                     "a" + a + "b" + b + "ga" + gamma + "di" + d)
            self.model = model
        # The case where parameters are controlled and there is an initial disturbance
        elif diff_para != None and controlled_inital == True:
            a, b, gamma, d = diff_para.convert_to_str()
            self.path_to_msh_file = "../../Gmsh/" + geometry + "/" + n_holes + ".msh"
            self.mesh_folder = "../../../Intermediate/" + geometry + "_mesh/" + n_holes + "_mesh/"
            self.pwd_folder = ("../../../Result/" + model + "/" + geometry + "/pwd_files/" + "h" + hole 
                              + "_d" + "r" + r + "x" + x_mid + "y" + y_mid +  "_k" 
                              + "a" + a + "b" + b + "ga" + gamma + "di" + d)
            self.file_save_folder = ("../../../Intermediate/" + model + "_files/" + geometry + "/h" + hole 
                                    + "_d" + "r" + r + "x" + x_mid + "y" + y_mid +  "_k" 
                                    + "a" + a + "b" + b + "ga" + gamma + "di" + d)
            self.model = model
        else:
            print("Error, improper file-locations entry")
            sys.exit(1)
        
    # Ensure that the save-locations are folders, and contain information
    # about changes to the standard parameters
    def change_standard_para(self, para_diff_str):
        self.pwd_folder += para_diff_str
        self.file_save_folder += para_diff_str
        

# -----------------------------------------------------------------------------------
# Start of functions 
# -----------------------------------------------------------------------------------

# Function that will check if a model's parameters fulfill the classical Turing
# conditions.
# Args:
#     param, a parameter object
#     model, the model
# Returns:
#     True if the conditions are fulfilled, else false
def check_if_param_fulfill_turing(param, model):
    if model == "Gierer":
        u = (param.a + 1) / param.b
        v = u**2
        fu = -param.b + 2 * u/v
        fv = -1/u**2
        gu = 2*u
        gv = -1
    elif model == "Schankenberg":
        u = param.a + param.b
        v = param.b / (param.a + param.b)**2
        fu = 2*u*v - 1
        fv = u**2
        gu = -2*u*v
        gv = -u**2
    
    # Check if the conditions are fulfilled 
    cond1 = fu + gv < 0
    cond2 = fu*gv - fv*gu > 0
    cond3 = param.d*fu + gv > 0
    cond4 = fu*gv - fv*gu < (param.d*fu +gv)**2/(4*param.d)
    if cond1 and cond2 and cond3 and cond4:
        return True
    else:
        return False


# Function that will convert a mesh on .msh format to xdmf-format, and extract the
# subdomains, the entire mesh, and lines for a 2d geometry. The result will be stored
# in a used-provided folder.
# Args:
#    file_locations, an object of type file locations (contains all file locations for a model)
# Returns:
#    void
def read_and_convert_mesh(file_locations):
    # Convert the mesh to xdmf format 
    msh = meshio.read(file_locations.path_to_msh_file)
    
    # Create directory if doesn't exist 
    if not os.path.isdir(file_locations.mesh_folder):
        os.makedirs(file_locations.mesh_folder)
    
    # Write all triangles
    meshio.write(file_locations.mesh_folder + "mesh.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]}))
    
    # Write the triangle subdomains 
    meshio.write(file_locations.mesh_folder + "subdomains.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]},
                                                              cell_data={"triangle": {"name_to_read": msh.cell_data["triangle"]["gmsh:physical"]}}))
    
    # Write the lines 
    meshio.write(file_locations.mesh_folder + "lines.xdmf", meshio.Mesh(points=msh.points, cells={"line": msh.cells["line"]},
                                                         cell_data={"line": {"name_to_read": msh.cell_data["line"]["gmsh:physical"]}}))


# Function that will formulate the nonlinear system (F) for the Schankenberg model
# with a specific dx-measure. More specific it contains the implicit backward
# Euler formulation 
# Args:
#    param, a parameter class object
#    u_1, u_2, the trial functions
#    v_1, v_2, the test functions
#    u_n1, u_n2, trial functions for previous step
#    dt_inv, the inverted time
#    dx, the space measure
# Returns:
#    The variational formula for backward Euler 
def formulate_FEM_equation(param, u_1, u_2, v_1, v_2, u_n1, u_n2, dt_inv, dx):
    dt_inv = Constant(dt_inv)
    a = Constant(param.a)
    b = Constant(param.b)
    d = Constant(param.d)
    gamma = Constant(param.gamma)
    
    F = dt_inv*(u_1 - u_n1)*v_1*dx - gamma*(a - u_1 + (u_1**2)*u_2)*v_1*dx \
        + inner(grad(u_1), grad(v_1))*dx + dt_inv*(u_2 - u_n2)*v_2*dx \
        - gamma*(b - (u_1**2)*u_2)*v_2*dx + d*inner(grad(u_2), grad(v_2))*dx
    
    return F


# Function that will formulate the nonlinear system (F) for the Schankenberg model
# with a specific dx-measure. This approaches uses the explicit Euler method when
# stepping in time (but note that it is an implicit and explicit mixture).
# Args:
#    param, a parameter class object
#    u_1, u_2, the trial functions
#    v_1, v_2, the test functions
#    u_n1, u_n2, trial functions for previous step
#    dt_inv, the inverted time
#    dx, the space measure
# Returns:
#    The variational formula for forward Euler 
def formulate_FEM_schnakenberg(param, u_1, u_2, v_1, v_2, u_n1, u_n2, dt_inv, dx):
    
    # Model parameters 
    dt_inv = Constant(dt_inv)
    a = Constant(param.a)
    b = Constant(param.b)
    d = Constant(param.d)
    gamma = Constant(param.gamma)
    
    F = dt_inv*(u_1 - u_n1)*v_1*dx - gamma*(a - u_n1 + (u_n1**2)*u_n2)*v_1*dx \
        + inner(grad(u_1), grad(v_1))*dx + dt_inv*(u_2 - u_n2)*v_2*dx \
        - gamma*(b - (u_n1**2)*u_n2)*v_2*dx + d*inner(grad(u_2), grad(v_2))*dx
    
    return F

# Function that will formulate the nonlinear system (F) for the Gierer-meinhardt
# model with a specific dx-measure. This approaches uses the explicit Euler method when
# stepping in time (but note that it is an implicit and explicit mixture).
# Args:
#    param, a parameter class object
#    u_1, u_2, the trial functions
#    v_1, v_2, the test functions
#    u_n1, u_n2, trial functions for previous step
#    dt_inv, the inverted time
#    dx, the space measure
# Returns:
#    The variational formula for forward Euler 
def formulate_FEM_gierer(param, u_1, u_2, v_1, v_2, u_n1, u_n2, dt_inv, dx):
    
    # Model parameters
    dt_inv = Constant(dt_inv)
    a = Constant(param.a)
    b = Constant(param.b)
    d = Constant(param.d)
    gamma = Constant(param.gamma)
    
    F = dt_inv*(u_1 - u_n1)*v_1*dx - gamma*(a - b * u_n1 + (u_n1**2) / u_n2)*v_1*dx \
        + inner(grad(u_1), grad(v_1))*dx + dt_inv*(u_2 - u_n2)*v_2*dx \
        - gamma*(u_n1 ** 2 - u_n2)*v_2*dx + d*inner(grad(u_2), grad(v_2))*dx
    
    return F


# FEM-solver for the backward Euler method, note that currently 10-iterations are used
# when solving the system of equations.
# Args:
#     V, the function space
#     F, the variational formulation
#     u, the test functions
#     u_n, the previous step solution
#     file_locations, a folder with the file locations for a current model 
#     dt, the step size
#     n_time_step, the number of time steps used when solving the problem 
# Returns:
#     void 
def solve_backward_euler(V, F, u, u_n, dt, file_locations, n_time_step):
    
    # Files for saving the end result
    file1 = file_locations.pwd_folder + "/u_1.pvd"
    file2 = file_locations.pwd_folder + "/u_2.pvd"
    vtkfile_u_1 = File(file1)
    vtkfile_u_2 = File(file2)
    
    # Solving the problem in time
    t = 0
    dU = Function(V)
    for n in range(n_time_step):
        if (n + 1) % 10 == 0:
            print("Time step {} of {}".format(n+1, n_time_step))
        
        # Update time-step 
        t += dt
        
        # Iterate 10-times when trying to solve the problem
        # In the future look into explicit stepping 
        for m in range(0, 10):
            A = assemble(derivative(F, u), keep_diagonal = True) # Define the derivative
            R = assemble(F) # Define the right hand side
            A.ident_zeros() # Fix "zero"-rows to get consistent system
            solve(A, dU.vector(), R) # Solver matrix system
            u.assign(u - dU) # Update solution one step    
        
        _u_1, _u_2, = u.split()
        # Save current time-step to file (only save every tenth step)
        if n % 10 == 0 or n == 1:
            vtkfile_u_1 << (_u_1, t)
            vtkfile_u_2 << (_u_2, t)
        u_n.assign(u)


# FEM-solver for the forward Euler system. Besides solving the system the
# maximal u1 concentration is saved at each iterations and written to
# file. 
# Args:
#     F, the variational formulation
#     u, the test functions
#     u_n, the previous step solution
#     file_locations, a class object with the file locations for the model 
#     dt, the step size
#     n_time_step, the number of time steps used when solving the problem 
# Returns:
#     u_1, u_2, the states at the final time (can be used for plotting)
def solve_forward_euler(V, F, u_n, file_locations, dt, n_time_step):
    # Files for saving the end result
    file1 = file_locations.pwd_folder + "/u_1.pvd"
    file2 = file_locations.pwd_folder + "/u_2.pvd"
    vtkfile_u_1 = File(file1)
    vtkfile_u_2 = File(file2)
    
    # Solving the problem in time
    t = 0
    U = Function(V)
    U.assign(u_n)
    
    # Write the first time-step
    _u_1, _u_2, = U.split()
    vtkfile_u_1 << (_u_1, t)
    vtkfile_u_2 << (_u_2, t)
    
    # Solving a linear system, all non-linear into last-vector 
    FA = lhs(F)
    Fb = rhs(F)
    A = assemble(FA, keep_diagonal = True)
    A.ident_zeros()
    
    # Vector for keeping track on maximum concentration 
    max_conc = np.zeros((n_time_step, 2))
    
    for n in tqdm(range(n_time_step)):
        
        # Update time-step
        t += dt
        
        # Solve linear variational problem for time step
        b = assemble(Fb)
        solve(A,  U.vector(), b)        
        
        _u_1, _u_2, = U.split()
        # Save current time-step to file (only save every tenth step)
        if n % 10 == 0 or n == 1:
            vtkfile_u_1 << (_u_1, t)
            vtkfile_u_2 << (_u_2, t)
        
        # Store the maximum value
        max_conc[n, 0] = t
        u1_vec = _u_1.vector().get_local()
        # Avoid problems with two molecules being stored 
        len_vec = int(len(u1_vec) / 2)
        u1_vec = u1_vec[np.tile([True, False], len_vec)]
        max_conc[n, 1] = np.max(u1_vec)
        
        # Sanity check solution
        min_u1 = np.min(u1_vec)
        min_u2 = np.min(_u_2.vector().get_local())
        if min_u1 < 0:
            print("Error, negative u1 concentration")
            sys.exit(1)
        elif min_u2 < 0:
            print("Error, negative u2 concentration")
            sys.exit(1)
        
        # Update previous solution 
        u_n.assign(U)
        
    
    # Write the maximum concentration to file 
    data_to_save = pd.DataFrame({"time": max_conc[0:, 0], "Max_conc": max_conc[0:, 1]})
    file_save = file_locations.file_save_folder + "max_conc.csv"
    # Create folder if not present
    if not os.path.isdir(file_locations.file_save_folder):
        os.makedirs(file_locations.file_save_folder)
    # If the file doesn't exist write a header, else append the result 
    if not os.path.isfile(file_save):
        data_to_save.to_csv(file_save)
    else:
        data_to_save.to_csv(file_save, header=False, mode='a')
    # Returning the states at the final time
    return _u_1, _u_2


# Function that keep track of which parameter in the large region don't have
# standard value, and thus returns a string with information about this
# that is added to the file-name.
# For Schankenberg standard is: a=0.2, b=2.0, gamma=10, d=100
# For Gierer standard is: a=0.5, b=2.0, gamma=20, d=50
# Args:
#    file_locations, object that contains information of the model
#    param, the parameters
def check_param_against_standard(file_locations, param):
    str_to_return = ""
    a, b, gamma, d = param.a, param.b, param.gamma, param.d
    if file_locations.model == "Schankenberg":
        if a != 0.2 or b != 2.0 or gamma != 10.0 or d != 100:
            str_to_return += "_K"
            if a != 0.2:
                str_to_return += "a" + str(a).replace(".", "D")
            if b != 2.0:
                str_to_return += "b" + str(b).replace(".", "D")
            if gamma != 10:
                str_to_return += "ga" + str(gamma).replace(".", "D")
            if d != 100:
                str_to_return += "d" + str(d).replace(".", "D")
    elif file_locations.model == "Gierer":
        if a != 0.5 or b != 2.0 or gamma != 20.0 or d != 50:
            str_to_return += "_K"
            if a != 0.5:
                str_to_return += "a" + str(a).replace(".", "D")
            if b != 2.0:
                str_to_return += "b" + str(b).replace(".", "D")
            if gamma != 20.0:
                str_to_return += "ga" + str(gamma).replace(".", "D")
            if d != 50:
                str_to_return += "di" + str(d).replace(".", "D")
    str_to_return += "/"
    return str_to_return


# Function that will solve the Schankenberg reaction diffusion system when the
# when the holes are defined via subdomains with zero flux and that the concentration
# is zero within the holes. The mesh, with subdomains, is read from msh-file that
# is converted to xdmf. Note that a set of standard parameters is assumed:
# For Schankenberg: a=0.2, b=2.0, gamma=10, d=100
# For Gierer: a=0.5, b=2.0, gamma=20, d=50
# If different parameters are supplied it will be noted in the file-save names
# Args:
#     param, a parameter class object containing the parameters of the model
#     t_end, the end time when solving the system
#     n_time_step, the number of time steps
#     dx_index_list, a list of which space surface measures to solve the PDE over
#     file_locations, a class object containing all the file locations for saving the result
#     ic_parameters, class object containing the parameters for initial conditions.
#     diff_para, a different parameter object, this is done to adapt the parameters in certain regions 
#     use_backward, if true the backward Euler method is used for solving the problem,
#         by default explicit Euler is desired.
#     seed, the seed to use when generating the random initial states
# Returns:
#     pwd-files, the function will save and output pwd-files for plotting
#     max_conc-file, a csv-file containing the maximum concentration at each time step
#     t_end-file, a csv-file containing the concentrations at the end time. Note that upon running the function
#         several times the function will append already existing csv-files. 
def solve_fem(param, t_end, n_time_step, dx_index_list, file_locations, ic_par, diff_para, use_backward=False, seed=123):
    
    # Setting the seed to reproduce the result 
    np.random.seed(seed)
    
    # Check if the Turing-conditions are fulfilled
    if not check_if_param_fulfill_turing(param, file_locations.model):
        print("The model does not fulfill Turing instability")
    
    # Reading the mesh into FeniCS 
    mesh = Mesh()
    with XDMFFile(file_locations.mesh_folder + "mesh.xdmf") as infile:
        infile.read(mesh)
    
    # Reading the triangle subdomains and storing them 
    sub_domains_pre = MeshValueCollection("size_t", mesh, 2)
    with XDMFFile(file_locations.mesh_folder + "subdomains.xdmf") as infile:
        infile.read(sub_domains_pre, "name_to_read")
    sub_domains = cpp.mesh.MeshFunctionSizet(mesh, sub_domains_pre)
    
    # Reading the line boundaries and storing them 
    line_domain_pre = MeshValueCollection("size_t", mesh, 1)
    with XDMFFile(file_locations.mesh_folder + "lines.xdmf") as infile:
        infile.read(line_domain_pre, "name_to_read")
    line_domain = cpp.mesh.MeshFunctionSizet(mesh, line_domain_pre)
    
    # Only display potential errors when solving the PDE-system 
    set_log_level(40)
    
    # List to hold the current fem-models
    fem_models = [formulate_FEM_schnakenberg, formulate_FEM_gierer]
    
    # Initial values for the model Schankenberg model
    if file_locations.model == "Schankenberg":
        u0_1 = param.a + param.b
        u0_2 = (u0_1 - param.a) / (u0_1**2)
        formulate_FEM = fem_models[0]
    elif file_locations.model == "Gierer":
        u0_1 = (param.a + 1) / param.b
        u0_2 = u0_1 ** 2
        formulate_FEM = fem_models[1]
    else:
        print("Error, invalid model name")
        sys.exit(1)
    
    # Fix the file-save locations, this will change file-names
    # if changes occurred to the standard parameters. 
    str_diff = check_param_against_standard(file_locations, param)
    file_locations.change_standard_para(str_diff)
    
    # Step size in t, use dt_inf for numerical precision 
    dt_inv = 1 / (t_end / n_time_step)
    dt = t_end / n_time_step
    
    # Standard Lagrange elements
    P1 = FiniteElement('P', mesh.ufl_cell(), 1)
    element = MixedElement([P1, P1])
    
    # Function space to approximate solution in 
    V = FunctionSpace(mesh, element)
    u_n = Function(V)
    
    # Expression for the initial conditions 
    u0_exp = ic(u0_1 = u0_1, u0_2 = u0_2, sigma = 0.05,
                init_param=ic_par, element = V.ufl_element())
    
    # Test functions, and initial values 
    v_1, v_2 = TestFunctions(V)
    u_n.interpolate(u0_exp)
    u_n1, u_n2 = split(u_n)
    
    # Try to get the coordinates
    n = V.dim()  
    d = mesh.geometry().dim()
    dof_coordinates = V.tabulate_dof_coordinates().reshape(n, d)
    dof_x = dof_coordinates[:, 0]                 
    dof_y = dof_coordinates[:, 1]
    
    # Define the different measures, note that:
    # dx -> surfaces (in this case)
    # ds -> exterior boundaries
    # dS -> interior boundaries
    dx = Measure("dx", domain=mesh, subdomain_data=sub_domains)
    ds = Measure("ds", domain=mesh, subdomain_data=line_domain)
    dS = Measure("dS", domain=mesh, subdomain_data=line_domain)
    
    # Solve the FEM-equations 
    F = 0
    if use_backward == True:
        print("Backward weak form")
        u = Function(V)
        u_1, u_2 = split(u)
        for i in dx_index_list:
            F += formulate_FEM_schnakenberg(param, u_1, u_2, v_1, v_2, u_n1, u_n2, dt_inv, dx(i))
        
        # Add the domains with different parameter-values
        if diff_para != None:
            for i in diff_para.dx_list:
                F += formulate_FEM(diff_para.param, u_1, u_2, v_1, v_2, u_n1, u_n2, dt_inv, dx(i))
        
        # Solving the backward system 
        solve_backward_euler(V, F, u, u_n, dt, folder_save, n_time_step)
        
    elif use_backward == False:
        u_1, u_2 = TrialFunctions(V)
        for i in dx_index_list:
            F += formulate_FEM(param, u_1, u_2, v_1, v_2, u_n1, u_n2, dt_inv, dx(i))
        
        # Add the domains with different parameter-values
        if diff_para != None:
            for i in diff_para.dx_list:
                F += formulate_FEM(diff_para.param, u_1, u_2, v_1, v_2, u_n1, u_n2, dt_inv, dx(i))
        
        # Solve the PDE-system 
        u1, u2 = solve_forward_euler(V, F, u_n, file_locations, dt, n_time_step)
    
    # Write the t_end result to file 
    file_save = file_locations.file_save_folder + "t_end_data.csv"
    # If file doesn't exist write header, else append file 
    if not os.path.isfile(file_save):
        t_end_data = pd.DataFrame({"x": dof_x, "y": dof_y, "u1": u1.vector().get_local(), "id" : 1})
        # Fix problem with index (know which vector is saved)
        n_rep = int(len(t_end_data.index) / 2); id_mol = np.tile([1, 2], n_rep)
        t_end_data["id_mol"] = id_mol
        t_end_data.to_csv(file_save)
    else:
        # Increment maximum id
        df = pd.read_csv(file_save)
        id_new = np.max(df.loc[:, "id"]) + 1
        t_end_data = pd.DataFrame({"x": dof_x, "y": dof_y, "u1": u1.vector().get_local(), "id" : id_new})
        # Fix problem with index (know which molecule is which)
        n_rep = int(len(t_end_data.index) / 2); id_mol = np.tile([1, 2], n_rep)
        t_end_data["id_mol"] = id_mol
        t_end_data.to_csv(file_save, header=False, mode='a')
    


# Function that will solve the FEM-problem using the explicit-implicit solver but only
# output the pwd-files. The main goal of this function is to sanity check the result
# by creating several pwd-files, where every tenth file will be saved.
# Args:
#    param, the model parameters (class object)
#    mesh_folder, the path to the mesh folder
#    t_opt, a class object containing the t-options
#    model, the model that is run
#    init_par, the initial value parameters 
#    pwd_folder, the folder where the pwd-result is stored
#    seed, a user set seed
def fem_sanity_check(param, mesh_folder, t_opt, ic_par, pwd_folder, model, seed, diff_para=None):
    
    # Setting the seed to reproduce the result 
    np.random.seed(seed)
    
    dx_index_list = [1]
    t_end = t_opt.t_end
    n_time_step = t_opt.n_time_step
    
    # Reading the mesh into FeniCS 
    mesh = Mesh()
    with XDMFFile(mesh_folder + "mesh.xdmf") as infile:
        infile.read(mesh)
    
    # Reading the triangle subdomains and storing them 
    sub_domains_pre = MeshValueCollection("size_t", mesh, 2)
    with XDMFFile(mesh_folder + "subdomains.xdmf") as infile:
        infile.read(sub_domains_pre, "name_to_read")
    sub_domains = cpp.mesh.MeshFunctionSizet(mesh, sub_domains_pre)
    
    # Reading the line boundaries and storing them 
    line_domain_pre = MeshValueCollection("size_t", mesh, 1)
    with XDMFFile(mesh_folder + "lines.xdmf") as infile:
        infile.read(line_domain_pre, "name_to_read")
    line_domain = cpp.mesh.MeshFunctionSizet(mesh, line_domain_pre)
    
    # Only display potential errors when solving the PDE-system 
    set_log_level(40)
    
    # List to hold the current fem-models
    fem_models = [formulate_FEM_schnakenberg, formulate_FEM_gierer]
    
    # Initial values for the model Schankenberg model
    if model == "Schankenberg":
        u0_1 = param.a + param.b
        u0_2 = (u0_1 - param.a) / (u0_1**2)
        formulate_FEM = fem_models[0]
    elif model == "Gierer":
        u0_1 = (param.a + 1) / param.b
        u0_2 = u0_1 ** 2
        formulate_FEM = fem_models[1]
    else:
        print("Error, invalid model name")
        sys.exit(1)
    
    # Step size in t, use dt_inf for numerical precision 
    dt_inv = 1 / (t_end / n_time_step)
    dt = t_end / n_time_step
    
    # Standard Lagrange elements
    P1 = FiniteElement('P', mesh.ufl_cell(), 1)
    element = MixedElement([P1, P1])
    
    # Function space to approximate solution in 
    V = FunctionSpace(mesh, element)
    u_n = Function(V)
    
    # Expression for the initial conditions 
    u0_exp = ic(u0_1 = u0_1, u0_2 = u0_2, sigma = 0.05,
                init_param=ic_par, element = V.ufl_element())
    
    # Test functions, and initial values 
    v_1, v_2 = TestFunctions(V)
    u_n.interpolate(u0_exp)
    u_n1, u_n2 = split(u_n)
    
    dx = Measure("dx", domain=mesh, subdomain_data=sub_domains)
    ds = Measure("ds", domain=mesh, subdomain_data=line_domain)
    dS = Measure("dS", domain=mesh, subdomain_data=line_domain)
    
    # Solve the FEM-equations 
    F = 0
    
    u_1, u_2 = TrialFunctions(V)
    for i in dx_index_list:
        F += formulate_FEM(param, u_1, u_2, v_1, v_2, u_n1, u_n2, dt_inv, dx(i))
    
    # Add sub-regions with different parameter values 
    if diff_para != None:
        for i in diff_para.dx_list:
            F += formulate_FEM(diff_para.param, u_1, u_2, v_1, v_2, u_n1, u_n2, dt_inv, dx(i))
    
    # Solve the actual problem
    file1 = pwd_folder + "/u_1.pvd"
    file2 = pwd_folder + "/u_2.pvd"
    vtkfile_u_1 = File(file1)
    vtkfile_u_2 = File(file2)
    
    # Solving the problem in time
    t = 0
    U = Function(V)
    U.assign(u_n)
    
    # Write the first time-step
    _u_1, _u_2, = U.split()
    vtkfile_u_1 << (_u_1, t)
    vtkfile_u_2 << (_u_2, t)
    
    # Solving a linear system, all non-linear into last-vector 
    FA = lhs(F)
    Fb = rhs(F)
    A = assemble(FA, keep_diagonal = True)
    A.ident_zeros()
    
    for n in tqdm(range(n_time_step)):
        
        # Update time-step
        t += dt
        
        # Solve linear variational problem for time step
        b = assemble(Fb)
        solve(A,  U.vector(), b)        
        
        _u_1, _u_2, = U.split()
        # Save current time-step to file (only save every tenth step)
        if n % 10 == 0 or n == 1:
            vtkfile_u_1 << (_u_1, t)
            vtkfile_u_2 << (_u_2, t)
        
        u_n.assign(U)
    

# Function that will run the find-mesh size experiment by solving the PDE on different
# size meshes.
# Args:
#    lc_list, a list of the lc-values, note 01 -> 0.1, 006 -> 0.06
#    n_holes_list, a list over the number of holes to use 
#    model, a string of which model is used
#    geometry, a string (circles or rectangles)
#    t_opt, an object that contain the time options 
#    param, parameters for the correct model 
def find_mesh_size(lc_list, n_holes_list, model, geometry, t_opt, param):
    # Fix the times
    n_time_step = t_opt.n_time_step
    t_end = t_opt.t_end
    
    # General arguments
    dx = [1]
    init_param = init_val_param(controlled=True, geom="Rectangles", r_circle=0.25)
    print("Model = {}, geom = {}".format(model, geometry))
    # Loop through the lc-list
    for hole in n_holes_list:
        print("Solving for {}".format(hole))
        for lc in lc_list:
            # Correct file locations
            file_loc = file_locations_class(hole, geometry, model, init_param, find_mesh_size=True, lc=lc)
            
            # Change the lc-value to the current one in the list
            # by creating a temporary file (and then replacing the old)
            path_geo_file = file_loc.path_to_msh_file[0:-3] + "geo"
            temp_file = path_geo_file[0:-4] + "_temp" + ".geo"
            with open(path_geo_file) as f:
                lines = f.readlines()
            lines[0] = 'lc = ' + lc + ';\n'
            with open(temp_file, 'w') as f:
                for item in lines:
                    f.write("%s" % item)
            os.remove(path_geo_file)
            os.rename(temp_file, path_geo_file)
            
            # Mesh the file
            path_gmsh = "~/gmsh/gmsh-4.4.1-Linux64/bin/gmsh"
            os.system(path_gmsh + " -2 " + path_geo_file + " 1> /dev/null") 
            
            # Convert mesh and solve 
            read_and_convert_mesh(file_loc)
            diff_para = None
            solve_fem(param, t_end, n_time_step, dx, file_loc, init_param, diff_para)


# Function that will solve the PDE:s for different number of holes for the
# provided geometry case. In order to compare models each model will be run with the
# same parameter set and time interval. Note that currently only two models are
# allowed, Gierer and the "Schankenberg" model.
# Args:
#    t_opt, a t_opt class object used to store time 
#    param, an object of parameter class
#    ic_par, the parameters for the initial value (can be a list for each number of holes)
#    geometry, a string of the geometry being solved
#    model, a string specifying which model to solve
#    hole_list, a list object containing the number of holes, not that ic_par most match this list 
#    ic_controlled, whatever or not the initial values are controlled (false by default)
#    diff_par, a different parameter object list, most match elements with hole_list (and if provided ic-list)
#    seed, the seed used for generating the different start-guesses. 
def solve_rd_system(t_opt, param, geometry="Rectangles", model="Schankenberg", hole_list=[],
                    ic_par="", ic_controlled=False, diff_par_list=None, seed=123):
    
    t_end = t_opt.t_end
    n_time_step = t_opt.n_time_step
    for i in range(len(hole_list)):
        
        # Fix the initial values 
        if isinstance(ic_par, list):
            ic_par_i = ic_par[i]
        else:
            ic_par_i = init_val_param()
            diff_par_i = None
        
        # Fix the different parameters
        if diff_par_list == None:
            diff_par_i = None
        else:
            diff_par_i = diff_par_list[i]
        
        print("Number of holes = {}".format(hole_list[i]))
        
        # Fix the mesh
        file_loc = file_locations_class(hole_list[i], geometry, model, ic_par_i, ic_controlled, diff_para=diff_par_i)
        read_and_convert_mesh(file_loc)
        
        # Solve the zero holes case 
        dx_index_list = [1]
        solve_fem(param, t_end, n_time_step, dx_index_list, file_loc, ic_par_i, diff_par_i, seed=seed)
    


# Function that will solve the PDE for a model several times
# Args:
#    param, a parameter class object (for the correct model)
#    t_opt, the time options (a class object)
#    model, the model to run
#    geom, the geometry of the object 
#    times_run, the number of times to run the analysis
#    hole_list, a list object enumerating the number of holes to use 
#    ic_list, a list with initial conditions 
def run_rd_sim(param, t_opt, model, geom, times_run, hole_list, ic_list=None, diff_par_list=None):
    
    # Create a seed list for making it reproducible
    np.random.seed(123)
    seed_list = np.random.randint(low=1, high=1000, size=times_run)
    
    # Solve the system, first case if initial aren't controlled
    if ic_list == None and diff_par_list == None:
        k = 1
        for seed in seed_list:
            print("Model = {}, k = {} / {}, geom = {}".format(model, k, times_run, geom))
            solve_rd_system(t_opt, param, geom, model, hole_list, seed=seed)
            k += 1
    # If the initial are controlled 
    elif ic_list != None and diff_par_list == None:
        k = 1
        for seed in seed_list:
            print("Model = {}, k = {} / {}, geom = {}".format(model, k, times_run, geom))
            solve_rd_system(t_opt, param, geom, model, hole_list, ic_par=ic_list, ic_controlled=True, seed=seed)
            k += 1
    # If different parameters are used 
    elif diff_par_list != None and ic_list == None:
        k = 1
        for seed in seed_list:
            print("Model = {}, k = {} / {}, geom = {}".format(model, k, times_run, geom))
            solve_rd_system(t_opt, param, geom, model, hole_list, diff_par_list=diff_par_list, seed=seed)
            k += 1
    # If different parameters and controlled disturbance occurs
    elif diff_par_list != None and ic_list != None:
        k = 1
        for seed in seed_list:
            print("Model = {}, k = {} / {}, geom = {}".format(model, k, times_run, geom))
            solve_rd_system(t_opt, param, geom, model, hole_list, ic_par=ic_list, diff_par_list=diff_par_list, ic_controlled=True, seed=seed)
            k += 1


# Function that will perform the sanity check of the long-time simulation results by using random seeds
# given by the time-command. Note that both models are checked. Note that param_list etc contains
# the parameters for both model (Schankenberg and Gierer), where Schankenberg is the first model, hence
# input argument order is important 
# Args:
#     param_list, the model parameters
#     t_opt_list, the time options
#     times_run, the number of times to run the sanity check
#     tag, a tag for keeping tab on the experiment 
#     ic_par, the initial condition parameters
#     par_diff_list, the different parameters (None by default)
#     geom, the geometry (default circle)
#     n_holes, the number of holes (default five)
def sanity_check(param_list, t_opt_list, times_run, ic_par, tag, par_diff_list=None, geom="Circles", n_holes="Five_holes"):
    
    # Set the file-locations
    model_list = ["Schankenberg", "Gierer"]
    j = 0
    for model in model_list:
        file_loc = file_loc_sanity_check(geom, n_holes, model, tag)
        param = param_list[j]
        t_opt = t_opt_list[j]
        if par_diff_list != None:
            diff_para = par_diff_list[j]
        else:
            diff_para = None
        j += 1 
        
        print("Model = {}, tag = {}".format(model, tag))
        for i in range(times_run):
            # Random time-seed
            seed = int(time.time())
            pwd_folder = file_loc.pwd_folder + str(int(i+1))
            fem_sanity_check(param, file_loc.mesh_folder, t_opt, ic_par, pwd_folder,
                             model=file_loc.model, seed=seed, diff_para=diff_para)

