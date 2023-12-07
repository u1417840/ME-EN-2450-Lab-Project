#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  2 13:32:59 2023

@author: lukasgoulet
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io

from Lab9_rungekutta import runge_kutta
from Lab9_slir import slir

def plot_slir(beta_max, mu_L_min):
    
    mat_data = scipy.io.loadmat('EnvironmentalForcing.mat')
    T = mat_data["T"][0]
    tspan = mat_data["tspan"][0]

    T_beta = np.zeros(len(T))
    #populates the T_beta vector using T
    for i in range(len(T)):
        if (T[i] <= 0 or T[i] >= 35):
            T_beta[i] = 0
        else:
            T_beta[i] = 0.000241 * T[i]**2.06737 * (35 - T[i])**0.72859 

    # beta_max = 1
    beta = beta_max * T_beta

    # mu_L_min = 6
    mu_L = np.zeros(len(tspan))
    j = 0
    #populates the mu_L vector using T_beta
    for i in range(len(tspan)):
        mu_L[i] = sum(T_beta[j:i+1])
        while mu_L[i] > mu_L_min:
            j += 1
            mu_L[i] = sum(T_beta[j:i+1])

    mu_L_inv = 1 / mu_L
    mu_I = 10
    mu_I_inv = 1 / mu_I
    e = 0.001
    A_p = 5000
    params = [beta, mu_L_inv, mu_I_inv, e, T, A_p]

    B_i = 1
    P_i = 1.33 * 30 * (-0.35968 + 0.10789 * 15 - 0.00214 * 15 * 15) * 30
    S_i = P_i / A_p
    L_i = 0.01 * S_i
    I_i = 0
    R_i = mu_I * I_i    
    y_vec = np.array([B_i, P_i, S_i, L_i, I_i, R_i])

    #solves for t and y
    odefun = lambda t, y: slir(t, y, tspan, params)
    
    (B, P, S, L, I, R) = runge_kutta(odefun, tspan, y_vec)

    #plots the slir model for each variable
    plt.plot(tspan, P/A_p, label='Total Population', color='lime', linestyle = '-')
    plt.plot(tspan, B/A_p, label='Berry Population', color='magenta', linestyle = '--')
    plt.plot(tspan, S, label='Susceptible', color='black', linestyle = '-.')
    plt.plot(tspan, L, label='Latent', color='cyan', linestyle = '-')
    plt.plot(tspan, I, label='Infected', color='blue', linestyle = ':')
    plt.plot(tspan, R, label='Removed', color='red', linestyle = '-.')
    plt.legend()
    plt.xlim(0, 61)
    plt.ylim(0, 1.6)
    plt.xlabel('time (days)')
    plt.ylabel('Population (fraction of initial)')
    plt.title(f'SLIRP Model for Plant Disease Simulation Plots \n beta_max = {beta_max}, mu_L_min = {mu_L_min}')
    plt.show()
    
    return