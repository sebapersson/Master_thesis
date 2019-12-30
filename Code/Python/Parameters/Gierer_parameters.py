#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os

'''
This file will create a plot of the Turing-space for the Gierer-Meinhardt reaction 
diffusion model.
'''

# Function that will calculate the upper and lower limit for the parameter a
# Args:
#    u0, the starting concentration of the u species in the model
#    d, the diffusion ratio
# Returns:
#    a_low and a_upper, the upper and higher limits for a
def calc_a_upper(u0, d):
    a_upper = np.minimum(1 +  (d*u0 - np.sqrt(8*d**3*u0))/d**2, 1-u0/d)
    return a_upper

# Function that will calculate the lower limit for the parameter a
# Args:
#    u0, the starting concentration of the u species in the model
#    d, the diffusion ratio
# Returns:
#    a_low and a_upper, the upper and higher limits for a
def calc_a_lower(u0, d):
    a_lower = 1 - u0
    return a_lower


# -------------------------------------------------------------------------
# Calculating the limits 
# -------------------------------------------------------------------------

# Limits given d = 25
d = 25
u0_25 = np.linspace(start=0.001, stop=4.1, num = 1000)
a_25 = calc_a_upper(u0_25, d)
b_25 = (a_25 + 1) / u0_25

# Calculate the lower values 
a_lower = calc_a_lower(u0_25, d)
b_lower = (a_lower + 1) / u0_25

# Limits given d = 50
d = 50
u0_50 = np.linspace(start=0.0001, stop=4.1, num = 1000)
a_50 = calc_a_upper(u0_50, d)
b_50 = (a_50 + 1) / u0_50

# Limits given d = 100
d = 100
u0_100 = np.linspace(start=0.0001, stop=5.1, num = 1000)
a_100 = calc_a_upper(u0_100, d)
b_100 = (a_100 + 1) / u0_100

# Plotting the different curves
plt.style.use(plt.style.available[17])
f = plt.figure(figsize=(9, 6))
col = '#000000'
plt.plot(b_lower, a_lower, linewidth=2.5, color = col)
plt.fill_between(b_lower, a_lower + 0.03, a_lower, alpha = 0.2, color = col)
plt.plot(b_lower, a_lower, linewidth=2.5, color = '#d73027', linestyle='--')
plt.plot(b_25, a_25, linewidth=2.5, color = col)
plt.fill_between(b_25, a_25 - 0.03, a_25, alpha = 0.2, color = col)
plt.plot(b_50, a_50, linewidth=2.5, color = col)
plt.fill_between(b_50, a_50 - 0.03, a_50, alpha = 0.2, color = col)
plt.plot(b_100, a_100, linewidth=2.5, color = col)
plt.fill_between(b_100, a_100 - 0.03, a_100, alpha = 0.2, color = col)

plt.ylim((0, 0.8))
plt.xlim((0.3, 4))
plt.xlabel("b", fontsize = 16)
plt.ylabel("a", fontsize = 16)
plt.text(1.5, 0.15, "Lower limit", fontsize = 13)
plt.text(1.5, 0.53, "d = 25", fontsize = 13)
plt.text(1.5, 0.65, "d = 50", fontsize = 13)
plt.text(1.5, 0.75, "d = 100", fontsize = 13)

# Save the figure to disk
path_dir = "../../../Result/Param_space/"
if not os.path.isdir(path_dir):
    os.mkdir(path_dir)

path_save = path_dir + "Gierer_turing_space.pdf"
f.savefig(path_save, bbox_inches='tight')

