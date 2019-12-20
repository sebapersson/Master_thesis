#!/usr/bin/env python
from Simulate_RD_models import *

# ------------------------------------------------------------------------------
# Check for optimal mesh-size 
# ------------------------------------------------------------------------------
# Find mesh size Schankenberg
print("Checking for optimal mesh-size")
lc_list = ["0.1", "0.06", "0.08"];
n_holes_list = ["Zero_holes", "Five_holes", "Twenty_holes"]
t_opt = t_opt_class(7.5, 1500)
param = param_schankenberg(gamma=10, d=100)
#find_mesh_size(lc_list, n_holes_list, "Schankenberg", "Rectangles", t_opt, param)
#find_mesh_size(lc_list, n_holes_list, "Schankenberg", "Circles", t_opt, param)

# Find the mesh size Gierer 
t_opt = t_opt_class(1.5, 2000)
param = param_gierer(b = 2.0, a = 0.5, gamma = 20, d = 50)
#find_mesh_size(lc_list, n_holes_list, "Gierer", "Rectangles", t_opt, param)
#find_mesh_size(lc_list, n_holes_list, "Gierer", "Circles", t_top, param)

# ------------------------------------------------------------------------------
# Running simulations without a specific disturbance in the steady state
# ------------------------------------------------------------------------------
print("Solving without controlled disturbance")
hl = ["Zero_holes", "Five_holes", "Twenty_holes"]
# Schankenberg
param = param_schankenberg(gamma=10, d=100)
t_opt = t_opt_class(7.5, 2000)
run_rd_sim(param, t_opt, model="Schankenberg", geom="Rectangles", hole_list=hl, times_run=20)
run_rd_sim(param, t_opt, model="Schankenberg", geom="Circles", hole_list=hl, times_run=20)

# Gierer
param = param_gierer(b = 2.0, a = 0.5, gamma = 20, d = 50)
t_opt = t_opt_class(1.5, 2000)
run_rd_sim(param, t_opt, model="Gierer", geom="Rectangles", hole_list=hl, times_run=20)
run_rd_sim(param, t_opt, model="Gierer", geom="Circles", hole_list=hl, times_run=20)

# ------------------------------------------------------------------------------
# Running simulations with a controlled disturbance steady state 
# ------------------------------------------------------------------------------
print("Running with controlled initial disturbance")
# The initial controlled conditions 
ic_par_zero = init_val_param(True, 0, 0, 1.0)
ic_par_five = init_val_param(True, -0.15, 0.5, 0.25)
ic_par_twenty = init_val_param(True, 0.55, 1.05, 0.25)
ic_list = [ic_par_zero, ic_par_five, ic_par_twenty]

# Schankenberg
param = param_schankenberg(gamma=10, d=100)
t_opt = t_opt_class(7.5, 2000)
run_rd_sim(param, t_opt, model="Schankenberg", geom="Rectangles", hole_list=hl, ic_list=ic_list, times_run=20)
run_rd_sim(param, t_opt, model="Schankenberg", geom="Circles", hole_list=hl, ic_list=ic_list, times_run=20)

# Gierer
param = param_gierer(b = 2.0, a = 0.5, gamma = 20, d = 50)
t_opt = t_opt_class(1.5, 2000)
run_rd_sim(param, t_opt, model="Gierer", geom="Rectangles", hole_list=hl, ic_list=ic_list, times_run=20)
run_rd_sim(param, t_opt, model="Gierer", geom="Circles", hole_list=hl, ic_list=ic_list, times_run=20)

# ------------------------------------------------------------------------------
# Running simulations with disturbed parameters 
# ------------------------------------------------------------------------------
hl = ["Five_holes", "Twenty_holes"]

# Schankenberg
param = param_schankenberg(gamma=10, d=100)
t_opt = t_opt_class(7.5, 2000)
diff_par = [diff_param_class([2], param_schankenberg(a=0.5, b=2.0, d=100, gamma=10)),
            diff_param_class([15], param_schankenberg(a=0.5, b=2.0, d=100, gamma=10))]
run_rd_sim(param, t_opt, model="Schankenberg", geom="Rectangles", hole_list=hl, diff_par_list=diff_par, times_run=20)

# Gierer
param = param_gierer(b = 2.0, a = 0.5, gamma = 20, d = 50)
t_opt = t_opt_class(1.5, 2000)
diff_par = [diff_param_class([2], param_gierer(a=0.6, b=2.0, d=50, gamma=20)),
            diff_param_class([15], param_gierer(a=0.6, b=2.0, d=50, gamma=20))]
run_rd_sim(param, t_opt, model="Gierer", geom="Rectangles", hole_list=hl, diff_par_list=diff_par, times_run=20)
