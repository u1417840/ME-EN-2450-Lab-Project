#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 15:56:56 2023

@author: lukasgoulet
"""

import matplotlib.pyplot as plt
import scipy.io

from Lab9_plots import plot_slir

mat_data = scipy.io.loadmat('EnvironmentalForcing.mat')
T = mat_data["T"][0]
tspan = mat_data["tspan"][0]

#part a
beta_max = 1
mu_L_min = 6
sol = plot_slir(beta_max, mu_L_min)

#part d
plt.plot(tspan, T)
plt.xlim(0, 61)
plt.xlabel('time (days)')
plt.ylabel('temperature (Celsius)')
plt.title('Forcing Air Temperature Data')
plt.show()

#part e
for beta_max in [0.5, 1, 5]:
    sol2 = plot_slir(beta_max, mu_L_min)
    
beta_max = 1
#part f
for mu_L_min in [3, 6, 20]:
    sol3 = plot_slir(beta_max, mu_L_min)


