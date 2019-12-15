#!/usr/bin/env python
from Simulate_RD_models import *

# This file aims to try to create a pattern at a certain location by creating an extra strong
# disturbance towards a new bud-site.

# The initial controlled conditions 
ic_par_zero = init_val_param(True, 0, 0, 1.0)
ic_par_five = init_val_param(True, -0.15, 0.5, 0.25)
ic_par_twenty = init_val_param(True, 0.55, 1.05, 0.25)
ic_list = [ic_par_zero, ic_par_five, ic_par_twenty]

# The seed list for each case
np.random.seed(123)
times_run = 10
seed_list = np.random.randint(low=1, high=1000, size=times_run)

# ------------------------------------------------------------------------
# Schankenberg-model
# ------------------------------------------------------------------------
param = param_schankenberg(gamma=10, d=100)
t_end = 7.5
n_time_step = 1500

# Solve for the different geometries 
geometry = "Rectangles"
for seed in seed_list:
    solve_rd_system(n_time_step, t_end, param, ic_list, geometry, seed=seed, ic_controlled=True)
geometry = "Circles"
for seed in seed_list:
    solve_rd_system(n_time_step, t_end, param, ic_list, geometry, seed=seed, ic_controlled=True)

# Solving the circle case
# Run the code with different seeds
np.random.seed(123)
seed_list = np.random.randint(low=1, high=1000, size=times_run)
for seed in seed_list:
    solve_rd_system(n_time_step, t_end, param, ic_list, geometry, seed=seed, ic_controlled=True)

# ------------------------------------------------------------------------
# Gierer-model 
# ------------------------------------------------------------------------
# General parameters 
param = param_schankenberg(gamma=10, d=100)
t_end = 2.5
n_time_step = 2000

# Run for the different geometries 
geometry = "Rectangles"
for seed in seed_list:
    solve_rd_system(n_time_step, t_end, param, ic_list, geometry, model="Gierer", seed=seed, ic_controlled=True)

geometry = "Circles"
for seed in seed_list:
    solve_rd_system(n_time_step, t_end, param, ic_list, geometry, model="Gierer", seed=seed, ic_controlled=True)
