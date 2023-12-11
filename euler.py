import numpy as np

# change to be of the form
# [y] = TimeInt(odefun,t,dt,y0,DepFlux_sum(cnt),vine(cnt).mu_L)

def euler(odefun, t, dt, y0, e, g):
        
    slope = odefun(t, y0, e, g)
    y = y0 + slope * dt
    
    return y