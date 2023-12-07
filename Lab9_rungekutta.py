#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 13:10:37 2023

@author: lukasgoulet
"""

import numpy as np

#part c
def runge_kutta(odefun, tspan, y0):
    
    eqs = len(y0)
    n = len(tspan)
    y = np.zeros([eqs, n])
    h = tspan[1] - tspan[0]
    for i in range(eqs):
        y[i, 0] = y0[i]
    for j in range(1, n):
        k1 = np.array(odefun(j-1, y[:,j-1]))
        k2 = np.array(odefun(j-1, y[:,j-1] + h*k1/2))
        k3 = np.array(odefun(j, y[:,j-1] + h*k2/2))
        k4 = np.array(odefun(j, y[:,j-1] + h*k3))
        y[:,j] = y[:,j-1] + h * (k1 + 2*k2 + 2*k3 + k4)/6  
            
    P_b = y[0,:]
    P = y[1,:] 
    S = y[2,:]
    L = y[3,:]
    I = y[4,:]
    R = y[5,:]
    return (P_b, P, S, L, I, R)

