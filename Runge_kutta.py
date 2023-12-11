import numpy as np

# change to be of the form
# [y] = TimeInt(odefun,t,dt,y0,DepFlux_sum(cnt),vine(cnt).mu_L)

def rk4(odefun, t, dt, y0, DepFlux_sum(cnt),vine(cnt).mu_L):
    tLen = len(t)
    yLen = len(x0)
    yList = np.zeros((tLen, yLen))
    yList[0, :] = y0
    for i in range(len(t) - 1):
        h = t[i+1] - t[i]
        k1 = np.array(odefun(t[i], yList[i, :]))
        k2 = np.array(odefun(t[i] + 0.5 * dt, yList[i, :] + 0.5 * k1 * dt))
        k3 = np.array(odefun(t[i] + 0.5 * dt, yList[i, :] + 0.5 * k2 * dt))
        k4 = np.array(odefun(t[i] + h, yList[i, :] + k3 * dt))
        slope = (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0
        yList[i+1, :] = yList[i, :] + slope * dt
    return yList

        






        
    
    
