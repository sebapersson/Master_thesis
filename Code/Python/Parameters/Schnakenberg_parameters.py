#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

'''
This file will create a plot of the Turing-space for the Schnakenberg reaction 
diffusion model.
'''

# Function that will calculate the upper and lower limit for the parameter a
# Args:
#    u0, the starting concentration of the u species in the model
#    d, the diffusion ratio
# Returns:
#    a_low and a_upper, the upper and higher limits for a 
def calc_a_upper(u0, d):
    a_upper = u0/2 * (1 - 2*u0 / np.sqrt(d) - u0**2/d)
    return a_upper


# Calculate the lower limit 
u0_lower = np.linspace(start=0.0, stop=1.2, num = 100)
a_low = u0_lower/2 * (1 - u0_lower**2)
b_low = u0_lower - a_low

# Calculate the limit when d = 25
u0_25 = np.linspace(start=0.0, stop=2.2, num = 100)
a_25 = calc_a_upper(u0_25, 25)
b_25 = u0_25 - a_25

# Calculate the limit when d = 50
u0_50 = np.linspace(start=0.0, stop=3.1, num = 100)
a_50 = calc_a_upper(u0_50, 50)
b_50 = u0_50 - a_50

# Calculate the limit when d = 100
u0_100 = np.linspace(start=0.0, stop=4.3, num = 100)
a_100 = calc_a_upper(u0_100, 100)
b_100 = u0_100 - a_100

# For the inf case 
u0_inf = np.linspace(start=0.0, stop=1.3, num = 100)
a_inf = 0.5 * u0_inf
b_inf = u0_inf - a_inf

# Plotting the different curves
plt.plot(b_low, a_low, linewidth=2.5, color = '#0868ac')
plt.fill_between(b_low, a_low, a_low + 0.03, alpha = 0.4, color = '#0868ac')
plt.plot(b_25, a_25, linewidth=2.5, color = '#a8ddb5')
plt.fill_between(b_25, a_25 - 0.03, a_25, alpha = 0.4, color = '#a8ddb5')
plt.plot(b_50, a_50, linewidth=2.5, color = '#7bccc4')
plt.fill_between(b_50, a_50 - 0.03, a_50, alpha = 0.4, color = '#7bccc4')
plt.plot(b_100, a_100, linewidth=2.5, color = '#4eb3d3')
plt.fill_between(b_100, a_100 - 0.03, a_100, alpha = 0.4, color = '#4eb3d3')
plt.plot(b_inf, a_inf, linewidth=2.5, color = '#0868ac')
plt.fill_between(b_inf, a_inf - 0.03, a_inf, alpha = 0.4, color = '#0868ac')
plt.ylim((0, 0.6))

plt.xlabel("a", fontsize = 16)
plt.ylabel("b", fontsize = 16)
plt.text(0.55, 0.04, "Lower limit", fontsize = 12)
plt.text(1.5, 0.2, "d = 25", fontsize = 12)
plt.text(2.1, 0.32, "d = 50", fontsize = 12)
plt.text(2.75, 0.43, "d = 100", fontsize = 12)
plt.text(0.60, 0.55, "d = infinity", fontsize = 12)

plt.show()



