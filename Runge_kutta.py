import numpy as np

# change to be of the form
# [y] = TimeInt(odefun,t,dt,y0,DepFlux_sum(cnt),vine(cnt).mu_L)

def rk4(odefun, t, dt, y0, DepFlux_sum(cnt),vine(cnt).mu_L):
    
    eqs = len(y0)
    n = len(t)
    y = np.zeros([eqs, n])
    h = t[1] - t[0]
    for i in range(eqs):
        y[i, 0] = y0[i]
    for j in range(n - 1):
        k1 = np.array(odefun(j, y[:,j], DepFlux_sum(cnt),vine(cnt).mu_L))
        k2 = np.array(odefun(j + 0.5, y[:,j] + 0.5 * k1 * dt, DepFlux_sum(cnt),vine(cnt).mu_L))
        k3 = np.array(odefun(j + 0.5, y[:,j] + 0.5 * k2 * dt, DepFlux_sum(cnt),vine(cnt).mu_L))
        k4 = np.array(odefun(j + dt, y[:,j] + k3 * dt, DepFlux_sum(cnt),vine(cnt).mu_L))
        slope = (k1 + 2*k2 + 2*k3 + k4)/6  
        y[:,j + 1] = y[:,j] + slope * dt
            
    P_b = y[0,:]
    P = y[1,:] 
    S = y[2,:]
    L = y[3,:]
    I = y[4,:]
    R = y[5,:]
    return (P_b, P, S, L, I, R)






        
    
    
