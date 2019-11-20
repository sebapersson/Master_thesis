#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from dolfin import *
import meshio
import os

# -----------------------------------------------------------------------------------
# Classes for solving the problem
# -----------------------------------------------------------------------------------

# Class to hold the parameters for the Schankenberg model 
class param_schankenberg:
    def __init__(self,a=0.2,b=2.0,gamma=100.0,d=50.0):
        self.a = a
        self.b = b
        self.gamma = gamma
        self.d = d

# Class for constructing the initial conditions for each state (python syntax)
# Note that the pattern formation is sensitive to the noise sigma. 
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

# -----------------------------------------------------------------------------------
# Start of functions 
# -----------------------------------------------------------------------------------

# Function that will convert a mesh on .msh format to xdmf-format, and extract the
# subdomains, the entire mesh, and lines for a 2d geometry. The result will be stored
# in a used-provided folder.
# Args:
#    path_to_msh_file, path to the msh-file assuming the directory the script is run form
#    folder_save, path to the folder where the result will be saved.
# Returns:
#    void
def read_and_convert_mesh(path_to_msh_file, folder_save):
    # Convert the mesh to xdmf format 
    msh = meshio.read(path_to_msh_file)
    
    # Create directory if doesn't exist 
    if not os.path.isdir(folder_save):
        os.mkdir(folder_save)
    
    # Write all triangles
    meshio.write(folder_save + "mesh.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]}))
    
    # Write the triangle subdomains 
    meshio.write(folder_save + "subdomains.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]},
                                                              cell_data={"triangle": {"name_to_read": msh.cell_data["triangle"]["gmsh:physical"]}}))
    
    # Write the lines 
    meshio.write(folder_save + "lines.xdmf", meshio.Mesh(points=msh.points, cells={"line": msh.cells["line"]},
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

# -----------------------------------------------------------------------------------
# End of functions 
# -----------------------------------------------------------------------------------


# Function that will solve the Schankenberg reaction diffusion system when the
# when the holes are defined via subdomains with zero flux and that the concentration
# is zero within the holes. The mesh, with subdomains, is read from msh-file that
# is converted to xdmf.
# Args:
#     param, a parameter class object containing the parameters of the model
#     t_end, the end time when solving the system
#     n_time_step, the number of time steps
#     mesh_folder, path to the folder where the mesh is located
#     df_index_list, a list of which measures to use when defining the sub-domains
#     folder_save, the folder to save the result in 
def solve_schankenberg_sub_domain_holes(param, t_end, n_time_step, mesh_folder, folder_save):
    
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
    
    # Only display potential errors when solving the problem  
    set_log_level(40)
    
    # Initial values for the Schankenberg model 
    u0_1 = param.a + param.b
    u0_2 = (u0_1 - param.a) / (u0_1**2)
    
    # Step size in t, use dt_inf for numerical precision 
    dt_inv = 1 / (t_end / n_time_step)
    dt = t_end / n_time_step
    
    # Standard Lagrange elements
    P1 = FiniteElement('P', mesh.ufl_cell(), 1)
    element = MixedElement([P1, P1])
    
    # Function space to approximate solution in 
    V = FunctionSpace(mesh, element)
    u = Function(V)
    u_n = Function(V)
    
    # Expression for the initial conditions 
    u0_exp = ic(u0_1 = u0_1, u0_2 = u0_2, sigma = 0.05, element = V.ufl_element())
    
    # Test and trial functions 
    v_1, v_2 = TestFunctions(V)
    u_1, u_2 = split(u)
    u_n.interpolate(u0_exp)
    u_n1, u_n2 = split(u_n)
    
    # Define the different measures, note that:
    # dx -> surfaces (in this case)
    # ds -> exterior boundaries
    # dS -> interior boundaries
    dx = Measure("dx", domain=mesh, subdomain_data=sub_domains)
    ds = Measure("ds", domain=mesh, subdomain_data=line_domain)
    dS = Measure("dS", domain=mesh, subdomain_data=line_domain)
    
    # The weak form, backward Euler in time, solve over all 
    F = 0
    for i in dx_index_list:
        F += formulate_FEM_equation(param, u_1, u_2, v_1, v_2, u_n1, u_n2, dt_inv, dx(i))
        
    # Files for saving the end result
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
        
        # Update time-step 
        t += dt
        
        # Iterate 10-times when trying to solve the problem
        # In the future look  into explicit stepping 
        for m in range(0, 10):
            A = assemble(derivative(F, u), keep_diagonal = True) # Define the derivative
            R = assemble(F) # Define the right hand side
            A.ident_zeros() # Fix "zero"-rows to get consistent system
            solve(A, dU.vector(), R) # Solver matrix system
            u.assign(u - dU) # Update solution one step    
        
        # Save current time-step to file 
        _u_1, _u_2, = u.split()
        vtkfile_u_1 << (_u_1, t)
        vtkfile_u_2 << (_u_2, t)
        u_n.assign(u)


# -----------------------------------------------------------------------------------------
# Rectangle with one hole
# -----------------------------------------------------------------------------------------
# Parameters 
param = param_schankenberg(gamma = 10, d=100)
t_end = 3
n_time_step = 200

# The index for the relevant surface measure 
dx_index_list = [1]

# Read the msh and store resulting files in Intermediate 
path_to_msh_file = "../../Gmsh/sub_dom.msh"
mesh_folder = "../../../Intermediate/Rectangle_one_hole/"
read_and_convert_mesh(path_to_msh_file, mesh_folder)

# Solve the system and store the result in test_sub_save
folder_save = "test_sub_save"
solve_schankenberg_sub_domain_holes(param, t_end, n_time_step, mesh_folder, folder_save)

