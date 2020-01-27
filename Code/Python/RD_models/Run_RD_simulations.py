#!/usr/bin/env python
from Simulate_RD_models import *


# ==================================================================================
# Check for optimal mesh-size
# ==================================================================================
# Find mesh size Schankenberg
print("Checking for optimal mesh-size")
lc_list = ["0.06", "0.03", "0.04"]
n_holes_list = ["Zero_holes", "Five_holes", "Twenty_holes"]
n_holes_list = ["Twenty_holes"]
t_opt = t_opt_class(7.5, 1500)
param = param_schankenberg(gamma=10, d=100)
find_mesh_size(lc_list, n_holes_list, "Schankenberg", "Rectangles", t_opt, param)
find_mesh_size(lc_list, n_holes_list, "Schankenberg", "Circles", t_opt, param)

# Find the mesh size Gierer 
t_opt = t_opt_class(1.5, 2000)
lc_list = ["0.06"]
param = param_gierer(b = 2.0, a = 0.5, gamma = 20, d = 50)
find_mesh_size(lc_list, n_holes_list, "Gierer", "Rectangles", t_opt, param)
find_mesh_size(lc_list, n_holes_list, "Gierer", "Circles", t_opt, param)

# ==================================================================================
# Running simulations without a specific disturbance in the steady state
# ==================================================================================
print("Running without controlled disturbance")
hl = ["Zero_holes", "Five_holes", "Twenty_holes"]
# Schankenberg
param = param_schankenberg(gamma=10, d=100)
t_opt = t_opt_class(7.5, 2000)
run_rd_sim(param, t_opt, model="Schankenberg", geom="Rectangles", hole_list=hl, times_run=1)
run_rd_sim(param, t_opt, model="Schankenberg", geom="Circles", hole_list=hl, times_run=5)

# Gierer
param = param_gierer(b = 2.0, a = 0.5, gamma = 20, d = 50)
t_opt = t_opt_class(1.5, 2000)
run_rd_sim(param, t_opt, model="Gierer", geom="Rectangles", hole_list=hl, times_run=5)
run_rd_sim(param, t_opt, model="Gierer", geom="Circles", hole_list=hl, times_run=5)


# ==================================================================================
# Running simulations with a controlled disturbance steady state
# ==================================================================================
print("Running with controlled initial disturbance")
# The initial controlled conditions 
ic_par_zero = init_val_param(True, "Circle", 0, 0, 0.25)
ic_par_five = init_val_param(True, "Circle", -0.15, 0.5, 0.25)
ic_par_twenty = init_val_param(True, "Circle", 0.55, 1.05, 0.25)
ic_list = [ic_par_zero, ic_par_five, ic_par_twenty]

# Schankenberg
param = param_schankenberg(gamma=10, d=100)
t_opt = t_opt_class(7.5, 2000)
run_rd_sim(param, t_opt, model="Schankenberg", geom="Rectangles", hole_list=hl, ic_list=ic_list, times_run=5)
run_rd_sim(param, t_opt, model="Schankenberg", geom="Circles", hole_list=hl, ic_list=ic_list, times_run=5)

# Gierer
param = param_gierer(b = 2.0, a = 0.5, gamma = 20, d = 50)
t_opt = t_opt_class(1.5, 2000)
run_rd_sim(param, t_opt, model="Gierer", geom="Rectangles", hole_list=hl, ic_list=ic_list, times_run=5)
run_rd_sim(param, t_opt, model="Gierer", geom="Circles", hole_list=hl, ic_list=ic_list, times_run=5)


# ==================================================================================
# Running simulations with disturbed parameters
# ==================================================================================
hl = ["Five_holes", "Twenty_holes"]
print("Running with different parameters in sub-region")

## Change parameters within the Turing-space, increase a 
# Schankenberg
param = param_schankenberg(gamma=10, d=100)
t_opt = t_opt_class(7.5, 2000)
# Increasing production 
diff_par = [diff_param_class([2], param_schankenberg(a=0.5, b=2.0, d=100, gamma=10)),
            diff_param_class([15], param_schankenberg(a=0.5, b=2.0, d=100, gamma=10))]
run_rd_sim(param, t_opt, model="Schankenberg", geom="Rectangles", hole_list=hl, diff_par_list=diff_par, times_run=5)
run_rd_sim(param, t_opt, model="Schankenberg", geom="Circles", hole_list=hl, diff_par_list=diff_par, times_run=5)

# Gierer
param = param_gierer(b = 2.0, a = 0.5, gamma = 20, d = 50)
t_opt = t_opt_class(1.5, 2000)
diff_par = [diff_param_class([2], param_gierer(a=0.6, b=2.0, d=50, gamma=20)),
            diff_param_class([15], param_gierer(a=0.6, b=2.0, d=50, gamma=20))]
run_rd_sim(param, t_opt, model="Gierer", geom="Rectangles", hole_list=hl, diff_par_list=diff_par, times_run=5)
run_rd_sim(param, t_opt, model="Gierer", geom="Circles", hole_list=hl, diff_par_list=diff_par, times_run=5)

# To save time, the following attempts tried only five holes
hl = ["Five_holes"]

## Increase diffusion (by a factor 2) and increase a-parameter
# Schankenberg
param = param_schankenberg(gamma=10, d=200)
t_opt = t_opt_class(7.5, 2000)
# Increasing production 
diff_par = [diff_param_class([2], param_schankenberg(a=0.6, b=2.0, d=200, gamma=10))]
run_rd_sim(param, t_opt, model="Schankenberg", geom="Rectangles", hole_list=hl, diff_par_list=diff_par, times_run=5)
run_rd_sim(param, t_opt, model="Schankenberg", geom="Circles", hole_list=hl, diff_par_list=diff_par, times_run=5)

# Gierer 
param = param_gierer(b = 2.0, a = 0.5, gamma = 20, d = 100)
t_opt = t_opt_class(1.5, 2000)
diff_par = [diff_param_class([2], param_gierer(a=0.65, b=2.0, d=100, gamma=20))]
run_rd_sim(param, t_opt, model="Gierer", geom="Rectangles", hole_list=hl, diff_par_list=diff_par, times_run=5)
run_rd_sim(param, t_opt, model="Gierer", geom="Circles", hole_list=hl, diff_par_list=diff_par, times_run=5)

## Heavily increase production (a parameter outside Turing-space)
# Schankenberg
param = param_schankenberg(gamma=10, d=100)
t_opt = t_opt_class(7.5, 2000)
diff_par = [diff_param_class([2], param_schankenberg(a=2.0, b=2.0, d=100, gamma=10))]
run_rd_sim(param, t_opt, model="Schankenberg", geom="Rectangles", hole_list=hl, diff_par_list=diff_par, times_run=5)
run_rd_sim(param, t_opt, model="Schankenberg", geom="Circles", hole_list=hl, diff_par_list=diff_par, times_run=5)

# Gierer 
param = param_gierer(b = 2.0, a = 0.5, gamma = 20, d = 50)
t_opt = t_opt_class(1.5, 2000)
diff_par = [diff_param_class([2], param_gierer(a=2.0, b=2.0, d=50, gamma=20))]
run_rd_sim(param, t_opt, model="Gierer", geom="Rectangles", hole_list=hl, diff_par_list=diff_par, times_run=5)
run_rd_sim(param, t_opt, model="Gierer", geom="Circles", hole_list=hl, diff_par_list=diff_par, times_run=5)

## Increase a and decrease b outside Turing-space
hl = ["Five_holes", "Seven_holes", "Twenty_holes"]
# Schankenberg, note original parameters are changed for param1 and param2
param_list = [param_schankenberg(gamma=5, d=800), param_schankenberg(gamma=10, d=800), param_schankenberg(gamma=10, d=100)]
t_opt = t_opt_class(1.0, 3000)
diff_par = [diff_param_class([2], param_schankenberg(a=2.5, b=0.2, d=100, gamma=10)),
            diff_param_class([2], param_schankenberg(a=2.5, b=0.2, d=100, gamma=10)),
            diff_param_class([15], param_schankenberg(a=2.5, b=0.2, d=100, gamma=10))]
for param in param_list:
    run_rd_sim(param, t_opt, model="Schankenberg", geom="Rectangles", hole_list=hl, diff_par_list=diff_par, times_run=5)
    run_rd_sim(param, t_opt, model="Schankenberg", geom="Circles", hole_list=hl, diff_par_list=diff_par, times_run=5)

# Gierer 
param = param_gierer(b = 2.0, a = 0.5, gamma = 20, d = 50)
t_opt = t_opt_class(1.5, 2000)
diff_par = [diff_param_class([2], param_gierer(a=2.0, b=0.5, d=50, gamma=20)),
            diff_param_class([2], param_gierer(a=2.0, b=0.5, d=50, gamma=20)),
            diff_param_class([15], param_gierer(a=2.0, b=0.5, d=50, gamma=20))]
run_rd_sim(param, t_opt, model="Gierer", geom="Rectangles", hole_list=hl, diff_par_list=diff_par, times_run=5)
run_rd_sim(param, t_opt, model="Gierer", geom="Circles", hole_list=hl, diff_par_list=diff_par, times_run=5)


# ==================================================================================
# Control both initial steady state and the parameters 
# ==================================================================================
hl = ["Five_holes", "Twenty_holes"]
print("Running with different parameters in sub-region and controlled initial")

# General parameters, disturbance at the same place 
ic_par_five = init_val_param(True, "Circle", 0.0, 0.0, 0.25)
ic_par_twenty = init_val_param(True, "Circle", 0.75, 0.85, 0.25)
ic_list = [ic_par_five, ic_par_twenty]

# Schankenberg, trying to see what happens with smaller holes 
param = param_schankenberg(gamma=20, d=50)
t_opt = t_opt_class(10.0, 2000)
diff_par = [diff_param_class([2], param_schankenberg(a=0.3, b=2.0, d=50, gamma=20)),
            diff_param_class([15], param_schankenberg(a=0.3, b=2.0, d=50, gamma=20))]
run_rd_sim(param, t_opt, model="Schankenberg", geom="Circles", hole_list=hl, ic_list=ic_list, diff_par_list=diff_par, times_run=5)
run_rd_sim(param, t_opt, model="Schankenberg", geom="Rectangles", hole_list=hl, ic_list=ic_list, diff_par_list=diff_par, times_run=5)

# Gierer
param = param_gierer(b = 2.0, a = 0.5, gamma = 20, d = 50)
t_opt = t_opt_class(1.5, 2000)
diff_par = [diff_param_class([2], param_gierer(a=0.6, b=2.0, d=50, gamma=20)),
            diff_param_class([15], param_gierer(a=0.6, b=2.0, d=50, gamma=20))]
run_rd_sim(param, t_opt, model="Gierer", geom="Rectangles", hole_list=hl, ic_list=ic_list, diff_par_list=diff_par, times_run=5)
run_rd_sim(param, t_opt, model="Gierer", geom="Circles", hole_list=hl, ic_list=ic_list, diff_par_list=diff_par, times_run=5)

# ==================================================================================
# Sanity check the solutions with a fine mesh 
# ==================================================================================
# General parameters
param_list = [param_schankenberg(gamma=10, d=100),
              param_gierer(b = 2.0, a = 0.5, gamma = 20, d = 50)]
t_opt_list = [t_opt_class(7.5, 2000), t_opt_class(1.5, 2000)]

# Disturbance in initial condition 
ic_par = init_val_param(True, "Circle", -0.15, 0.5, 0.25) # For five holes 
sanity_check(param_list, t_opt_list, 1, ic_par, tag="d_disturbed")

# Parameters disturbed in sub-region
ic_par = init_val_param()
diff_par_list = [diff_param_class([2], param_schankenberg(a=0.5, b=2.0, d=100, gamma=10)),
                 diff_param_class([2], param_gierer(a=0.6, b=2.0, d=50, gamma=20))]
sanity_check(param_list, t_opt_list, 4, ic_par, tag="k_disturbed", par_diff_list=diff_par_list)
