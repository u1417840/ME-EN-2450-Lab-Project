#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 00:19:27 2023

@author: lukasgoulet
"""

import numpy as np

def rk4(odefun, tspan, x0):
    tLen = len(tspan)
    xLen = len(x0)
    xList = np.zeros((tLen, xLen))
    xList[0, :] = x0
    for i in range(len(tspan) - 1):
        h = tspan[i+1] - tspan[i]
        k1 = np.array(odefun(tspan[i], xList[i, :]))
        k2 = np.array(odefun(tspan[i] + 0.5 * h, xList[i, :] + 0.5 * k1 * h))
        k3 = np.array(odefun(tspan[i] + 0.5 * h, xList[i, :] + 0.5 * k2 * h))
        k4 = np.array(odefun(tspan[i] + h, xList[i, :] + k3 * h))
        slope = (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0
        xList[i+1, :] = xList[i, :] + slope * h
    return tspan, xList

        






        
    
    
