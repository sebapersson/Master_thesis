#!/usr/bin/env python

import numpy as np
import pandas as pd
from dolfin import *
import meshio
import os
import sys 
from tqdm import tqdm

# -----------------------------------------------------------------------------------
# Classes for solving the problem
# -----------------------------------------------------------------------------------

# Class for holding the initial value parameters, if geom is equal to rectangle
# a quadratic rectangle will be disturbed with middle x_mid, y_mid and side length r_circle*2
class init_val_param:
    def __init__(self, controlled=False, geom="Circle", x_mid=0, y_mid=0, r_circle=0.25):
        self.controlled = controlled
        self.x_mid = x_mid
        self.y_mid = y_mid
        self.r_circle = r_circle
        self.geom = geom

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

# Class to hold the file-locations for a model, 
class file_locations_class:
    def __init__(self, n_holes, geometry, model, controlled_inital=False, find_mesh_size=False, lc="0"):
        # If the initial values aren't controlled
        if controlled_inital == False and find_mesh_size == False:
            self.path_to_msh_file = "../../Gmsh/" + geometry + "/" + n_holes + ".msh"
            self.mesh_folder = "../../../Intermediate/" + geometry + "_mesh/" + n_holes + "_mesh/"
            self.pwd_folder = "../../../Result/" + model + "/" + geometry + "_pwd_files/" + n_holes + "/"
            self.file_save_folder = "../../../Intermediate/" + model + "_files/" + geometry + "/" + n_holes + "_files/"
            self.model = model
        # Rename save-folders if initial values are controlled 
        elif controlled_inital == True and find_mesh_size == False:
            self.path_to_msh_file = "../../Gmsh/" + geometry + "/" + n_holes + ".msh"
            self.mesh_folder = "../../../Intermediate/" + geometry + "_mesh/" + n_holes + "_mesh/"
            self.pwd_folder = "../../../Result/" + model + "/" + geometry + "_pwd_files/" + n_holes + "_init_controlled/"
            self.file_save_folder = "../../../Intermediate/" + model + "_files/" + geometry + "_init_controlled/" + n_holes + "_files/"
            self.model = model
        # Rename save-folders if the aim is to find the mesh-size
        elif find_mesh_size == True:
            self.path_to_msh_file = "../../Gmsh/" + geometry + "/Rec_lc" + lc + ".msh"
            self.mesh_folder = "../../../Intermediate/" + geometry + "_mesh/" + "Find_mesh_size_lc" + lc + "_mesh/"
            self.pwd_folder = "../../../Result/" + model + "/" + geometry + "_pwd_files/" + "Find_mesh_size_lc_" + lc + "/"
            self.file_save_folder = "../../../Intermediate/" + model + "_files/" + geometry + "_find_mesh_size_lc/" + lc + "_files/"
            self.model = model
        else:
            print("Error, improper file-locations entry")
            sys.exit(1)

# -----------------------------------------------------------------------------------
# Start of functions 
# -----------------------------------------------------------------------------------

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
        max_conc[n, 1] = np.max(_u_1.vector().get_local())
        
        # Sanity check solution
        min_u1 = np.min(_u_1.vector().get_local())
        min_u2 = np.min(_u_2.vector().get_local())
        if min_u1 < 0 or min_u2 < 0:
            print("Error, negative concentration")
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


# Function that will solve the Schankenberg reaction diffusion system when the
# when the holes are defined via subdomains with zero flux and that the concentration
# is zero within the holes. The mesh, with subdomains, is read from msh-file that
# is converted to xdmf.
# Args:
#     param, a parameter class object containing the parameters of the model
#     t_end, the end time when solving the system
#     n_time_step, the number of time steps
#     dx_index_list, a list of which space surface measures to solve the PDE over
#     file_locations, a class object containing all the file locations for saving the result
#     ic_parameters, class object containing the parameters for initial conditions.
#     use_backward, if true the backward Euler method is used for solving the problem,
#         by default explicit Euler is desired.
#     seed, the seed to use when generating the random initial states
# Returns:
#     pwd-files, the function will save and output pwd-files for plotting
#     max_conc-file, a csv-file containing the maximum concentration at each time step
#     t_end-file, a csv-file containing the concentrations at the end time. Note that upon running the function
#         several times the function will append already existing csv-files. 
def solve_schankenberg_sub_domain_holes(param, t_end, n_time_step, dx_index_list, file_locations, ic_par,
                                        use_backward=False, seed=123):
    
    # Setting the seed to reproduce the result 
    np.random.seed(seed)
    
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
        
        # Solving the backward system 
        solve_backward_euler(V, F, u, u_n, dt, folder_save, n_time_step)
        
    elif use_backward == False:
        print("Solving using the forward Euler method")
        u_1, u_2 = TrialFunctions(V)
        for i in dx_index_list:
            F += formulate_FEM(param, u_1, u_2, v_1, v_2, u_n1, u_n2, dt_inv, dx(i))
        
        # Solve the PDE-system 
        u1, u2 = solve_forward_euler(V, F, u_n, file_locations, dt, n_time_step)
    
    # Write the t_end result to file 
    t_end_data = pd.DataFrame({"x": dof_x, "y": dof_y, "u1": u1.vector().get_local(), "u2": u2.vector().get_local()})
    file_save = file_locations.file_save_folder + "t_end_data.csv"
    # If file doesn't exist write header, else append file 
    if not os.path.isfile(file_save):
        t_end_data.to_csv(file_save)
    else:
        t_end_data.to_csv(file_save, header=False, mode='a')
    


# Function that will run the find-mesh size experiment by solving the PDE on different
# size meshes.
# Args:
#    lc_list, a list of the lc-values, note 01 -> 0.1, 006 -> 0.06
#    model, a string of which model is used
#    geometry, a string (circles or rectangles)
#    n_time_step, the number of time-steps
#    t_end, the end time
#    param, parameters for the correct model 
def find_mesh_size(lc_list, model, geometry, n_time_step, t_end, param):
    # General arguments
    dx = [1]
    init_param = init_val_param(controlled=True, geom=geometry, r_circle=0.25)
    # Loop through the lc-list
    for lc in lc_list:
        file_loc = file_locations_class("Zero", geometry, model, find_mesh_size=True, lc=lc)
        read_and_convert_mesh(file_loc)
        solve_schankenberg_sub_domain_holes(param, t_end, n_time_step, dx, file_loc, init_param)


# Function that will solve the PDE:s for different number of holes for the
# provided geometry case. In order to compare models each model will be run with the
# same parameter set and time interval. Note that currently only two models are
# allowed, Gierer and the "Schankenberg" model.
# Args:
#    n_time_step, the number of time-steps when solving the PDE
#    t_end, the end time when solving the pde
#    param, an object of parameter class
#    ic_par, the parameters for the initial value (can be a list for each number of holes)
#    geometry, a string of the geometry being solved
#    model, a string specifying which model to solve
#    ic_controlled, whatever or not the initial values are controlled (false by default)
#    seed, the seed used for generating the different start-guesses. 
def solve_rd_system(n_time_step, t_end, param, ic_par=init_val_param(), geometry="Rectangles", model="Schankenberg", ic_controlled=False, seed=123):
    
    # The file locations for each case
    file_locations_zero = file_locations_class("Zero_holes", geometry, model, ic_controlled)
    file_locations_five = file_locations_class("Five_holes", geometry, model, ic_controlled)
    file_locations_twenty = file_locations_class("Twenty_holes", geometry, model, ic_controlled)
    
    # Create all the different mesh
    read_and_convert_mesh(file_locations_zero)
    read_and_convert_mesh(file_locations_five)
    read_and_convert_mesh(file_locations_twenty)
    
    # Adapt unique initial conditions if ic_controlled is a list
    if isinstance(ic_par, list):
        ic_par_zero = ic_par[0]
        ic_par_five = ic_par[1]
        ic_par_twenty = ic_par[2]
    else:
        ic_par_zero = ic_par
        ic_par_five = ic_par
        ic_par_twenty = ic_par
    
    # Solve the zero holes case 
    print("Solving PDE " + geometry + " with zero holes")
    dx_index_list = [1]
    solve_schankenberg_sub_domain_holes(param, t_end, n_time_step, dx_index_list,
                                        file_locations_zero, ic_par_zero, seed=seed)
    
    # Solve the five holes case
    print("Solving PDE " + geometry + " with five holes")
    dx_index_list = [1]
    solve_schankenberg_sub_domain_holes(param, t_end, n_time_step, dx_index_list,
                                        file_locations_five, ic_par_five, seed=seed)
    
    # Solve the 20 holes case
    print("Solving PDE " + geometry + " with twenty holes")
    dx_index_list = [1]
    solve_schankenberg_sub_domain_holes(param, t_end, n_time_step, dx_index_list,
                                        file_locations_twenty, ic_par_twenty, seed=seed)


# Run the actual analysis 
def main():
    # -----------------------------------------------------------------------------------
    # Rectangle case
    # -----------------------------------------------------------------------------------
    # Solving the rectangle case
    param = param_schankenberg(gamma=10, d=100)
    t_end = 5
    n_time_step = 1000
    times_run = 20
    geometry = "Rectangles"
    # Run the code with different seeds
    np.random.seed(123)
    seed_list = np.random.randint(low=1, high=1000, size=times_run)
    for seed in seed_list:
        solve_rd_system(n_time_step, t_end, param, geometry, seed=seed)
    
    # -----------------------------------------------------------------------------------
    # Circle case 
    # -----------------------------------------------------------------------------------
    # Solving the circle case
    param = param_schankenberg(gamma=10, d=100)
    t_end = 5
    n_time_step = 1000
    times_run = 20
    geometry = "Circles"
    # Run the code with different seeds 
    np.random.seed(123)
    seed_list = np.random.randint(low=1, high=1000, size=times_run)
    for seed in seed_list:
        solve_rd_system(n_time_step, t_end, param, geometry, seed=seed)
    
    
    # -----------------------------------------------------------------------------------
    # Gierer 
    # -----------------------------------------------------------------------------------
    print("Gierer model")
    # Rectangle case 
    param = param_gierer(b = 2.0, a = 0.5, gamma = 20, d = 50)
    t_end = 1.5
    n_time_step = 2000
    geometry = "Rectangles"
    times_run = 20
    # Run the code with different seeds 
    np.random.seed(123)
    seed_list = np.random.randint(low=1, high=1000, size=times_run)
    for seed in seed_list:
        solve_rd_system(n_time_step, t_end, param, geometry, "Gierer", seed=seed)
    
    
    # Circle case 
    param = param_gierer(b = 2.0, a = 0.5, gamma = 20, d = 50)
    t_end = 1.5
    n_time_step = 2000
    geometry = "Circles"
    times_run = 20
    # Run the code with different seeds 
    np.random.seed(123)
    seed_list = np.random.randint(low=1, high=1000, size=times_run)
    for seed in seed_list:
        solve_rd_system(n_time_step, t_end, param, geometry, "Gierer", seed=seed)
    
    return 0


# Run the cases with controlling the mesh size
# Schankenberg model 
lc_list = ["01", "008", "006"]; geometry = "Rectangles"
n_time_step = 1500; t_end = 5.0; param = param_schankenberg(gamma=10, d=100)
#find_mesh_size(lc_list, "Schankenberg", geometry, n_time_step, t_end, param)

# Gierer model 
lc_list = ["01", "008", "006"]; geometry = "Rectangles"
n_time_step = 2000; t_end = 2.5
param = param_gierer(b = 2.0, a = 0.5, gamma = 20, d = 50)
find_mesh_size(lc_list, "Gierer", geometry, n_time_step, t_end, param)

