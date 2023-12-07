#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 18:41:58 2023

@author: lukasgoulet
"""

#part b
def slir(t, y, tspan, params):
    
    beta = params[0][t]
    mu_L_inv = params[1][t]
    mu_I_inv = params[2]
    e = params[3]
    T = params[4][t]
    A_p = params[5]
    
    B = y[0]
    S = y[2]
    L = y[3]
    I = y[4]
    P_b = B  #P_b = B_i
    tspan = tspan[t]

    T_E = -0.35968 + 0.10789 * T - 0.00214 * (T**2)
    dP_bdt = (0.1724 * P_b - 0.0000212 * P_b**2) * T_E
    tday = tspan + 30
    dP_ldt = (1.33 * tday) * T_E
    dPdt = dP_bdt + dP_ldt
    
    dSdt = (-beta * S * I) + (dPdt / A_p)
    dLdt = (beta * S * I) - (mu_L_inv * L) + e
    dIdt = (mu_L_inv * L) - (mu_I_inv * I)
    dRdt = (mu_I_inv * I)
    
    return (dP_bdt, dPdt, dSdt, dLdt, dIdt, dRdt)
