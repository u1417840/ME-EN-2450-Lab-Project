import numpy as np
from scipy.integrate import trapz
from Sall_temp_effect import sall_temp_effect

def latentperiod(istart, dt, Nsteps, mu_L_target, mu_L, T):
    # Latent period length as a function of Temperature (Calonnec et al. 2008)

    # Calculate temp effect array. This could be calculated once external to
    # the time loop for efficiency and passed to this function
    PT = np.zeros_like(T)
    for i in range(Nsteps):
        PT[i] = sall_temp_effect(T[i])  # Assuming you have a function sall_temp_effect

    # Now calculate the time to latent from istart (when a new vine was infected)
    flag = True
    ispan = 0

    for i in range(istart, Nsteps):
        if flag:
            PTint = trapz(PT[istart:i]) * dt
            if PTint >= mu_L_target:
                mu_L[istart:i] = dt * (i - istart)
                ispan = i
                istart = istart + 1
                flag = False
        else:
            PTint = trapz(PT[istart:i]) * dt
            if PTint >= mu_L_target:
                mu_L[ispan:i] = dt * (i - istart)
                istart = istart + 1
                ispan = i

    mu_L = 1.0 / mu_L
    # Set any infinite mu_L's back to zero
    infInd = np.isinf(mu_L)
    mu_L[infInd] = 0


    return mu_L
