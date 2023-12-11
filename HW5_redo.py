#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 00:19:27 2023

@author: lukasgoulet
"""

import numpy as np
import matplotlib.pyplot as plt

def ode(t, x):
    y, z = x
    return (-2 * y + 4 * np.exp(-t), (-y * (z ** 2)) / 3)

def euler(odefun, tspan, x0):
    tLen = len(tspan)
    xLen = len(x0)
    xList = np.zeros((tLen, xLen))
    xList[0, :] = x0
    for i in range(len(tspan) - 1):
        h = tspan[i+1] - tspan[i]
        slope = np.array(odefun(tspan[i], xList[i, :]))
        xList[i+1, :] = xList[i, :] + (slope * h)
    return tspan, xList

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
    
y0 = 2.0
z0 = 4.0
x0 = np.array([y0, z0])
t0 = 0
tf = 1
odefun = lambda t, x: ode(t, x)

h = 0.1
n = int((tf - t0) / h + 1)
tspan = np.linspace(t0, tf, n)
tList1, xList1 = euler(odefun, tspan, x0)
plt.plot(tList1, xList1[:,0], 'rx', label = f'eulers method with h = {h}')

tList2, xList2 = rk4(odefun, tspan, x0)
plt.plot(tList2, xList2[:,0], 'b.', label = f'eulers method with h = {h}')

h = 0.01
n = int((tf - t0) / h + 1)
tspan = np.linspace(t0, tf, n)
tList3, xList3 = euler(odefun, tspan, x0)
plt.plot(tList3, xList3[:,0], 'g--', label = f'eulers method with h = {h}')
plt.legend()
plt.show()


h = 0.1
n = int((tf - t0) / h + 1)
tspan = np.linspace(t0, tf, n)
tList4, xList4 = euler(odefun, tspan, x0)
plt.plot(tList4, xList4[:,1], 'rx', label = f'eulers method with h = {h}')

tList5, xList5 = rk4(odefun, tspan, x0)
plt.plot(tList5, xList5[:,1], 'b.', label = f'eulers method with h = {h}')

h = 0.01
n = int((tf - t0) / h + 1)
tspan = np.linspace(t0, tf, n)
tList6, xList6 = euler(odefun, tspan, x0)
plt.plot(tList6, xList6[:,1], 'g--', label = f'eulers method with h = {h}')
plt.legend()
plt.show()

def mass_spring(t, y, m, c, k):
    x, v = y
    return(v, (-k * x - c * v) / m)

def rk4(odefun, tspan, y0):
    tLen = len(tspan)
    tList = tspan
    yLen = len(y0)
    yList = np.zeros((tLen, yLen))
    yList[0, :] = y0
    for i in range(len(tList) - 1):
        h = tList[i + 1] - tList[i]
        k1 = np.array(odefun(tList[i], yList[i, :]))
        k2 = np.array(odefun(tList[i] + 0.5 * h, yList[i, :] + 0.5 * k1 * h))
        k3 = np.array(odefun(tList[i] + 0.5 * h, yList[i, :] + 0.5 * k2 * h))
        k4 = np.array(odefun(tList[i] + h, yList[i, :] + k3 * h))
        slope = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        yList[i + 1, :] = yList[i, :] + slope * h
    return tList, yList

m = 20
k = 20
cList = (5, 40, 200)
x0 = 1
v0 = 0 
y0 = ([x0, v0])
t0 = 0
tf = 15
n = int((tf - t0) / h + 1)
tspan = np.linspace(t0, tf, n)

odefun = lambda t, y: mass_spring(t, y, m, c, k)

for c in cList:
    tList, yList = rk4(odefun, tspan, y0)
    plt.plot(tList, yList[:, 0], label = f'rk4 with c = {c}')
    plt.legend()
        






        
    
    