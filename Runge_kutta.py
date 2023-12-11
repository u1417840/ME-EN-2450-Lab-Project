import numpy as np

def rk4(odefun, tspan, y0):
    tLen = len(tspan)
    yLen = len(x0)
    yList = np.zeros((tLen, yLen))
    yList[0, :] = y0
    for i in range(len(tspan) - 1):
        h = tspan[i+1] - tspan[i]
        k1 = np.array(odefun(tspan[i], yList[i, :]))
        k2 = np.array(odefun(tspan[i] + 0.5 * h, yList[i, :] + 0.5 * k1 * h))
        k3 = np.array(odefun(tspan[i] + 0.5 * h, yList[i, :] + 0.5 * k2 * h))
        k4 = np.array(odefun(tspan[i] + h, yList[i, :] + k3 * h))
        slope = (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0
        yList[i+1, :] = yList[i, :] + slope * h
    return yList

        






        
    
    
