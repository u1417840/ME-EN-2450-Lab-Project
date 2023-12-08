# SLIRPE function 
# % SLIRPE_model (pronounced like the 7-11 slushy drink) combines all the 
# % SLIR equations and the plant growth equations with an external source 
# % equation in one system. The dependent variables are stored in y as:
# % y(0) = B (amount of population, surface area, that is berries)
# % y(1) = P (total population, surface area, including berries and leaves: P = S+L+I+R)
# % y(2) = S (susceptible population)
# % y(3) = L (latent population)
# % y(4) = I (infectious population)
# % y(5) = R (recovered/removed population)
# % y(6) = E (amount of new infections from external sources)
# % y(7) = F (size of the spreading population, e.g. sporulating for a fungus) 
# %
# % and the parameters in a cell array:
# % p{0} = beta_max (max rate of colony growth/new infections)
# % p{1} = mu_i (inverse length of the infectious period in days)
# % p{2} = T (array of temperature in C)
# % p{3} = day (array of times in units of days)
# % p{4} = A (total plant surface area at reference time)
# % p{5} = Windspd (windspeed)
# % p{6} = Winddir (wind direction, currently not used here)
# % p{7} = eta     (release fraction scale factor)
# % p{8} = kappa   (release fraction scale factor)  
# % p{9}= xi      (release fraction offset)
# % p{10}= Gamma   (spore production multiple)
# % p{11}= alpha   (spore production 2nd factor)
# %
# % and the time input (idx) should be an integer for the iteration number
# % Note that E is not calculated by the function (only integrated in time)

import numpy as np
from Sall_temp_effect import sall_temp_effect

def SLIRPE_model(idx, y, e, mu_L, p):
    # Assign parameters
    beta_max, mu_I, T, day, A, Windspd, Winddir, eta, kappa, xi, Gamma, alpha = p[:12]

    # Assign variables
    B, P, S, L, I, R, E, F = y[:8]
    
    
    # Calculated parameters
    if np.ceil(idx) == np.floor(idx): #when we are at an interger step
        idx = int(idx)
        T_used = T[idx]
        day_used = day[idx]
        mu_L_used = mu_L[idx]
        m_used = Windspd[idx]
    else: #when we are at a half step (for rk4)
        idx = int(idx)
        T_used = 0.5*(T[idx]+T[idx+1])
        day_used = 0.5*(day[idx]+day[idx+1])
        mu_L_used = 0.5*(mu_L[idx]+mu_L[idx+1])
        m_used = 0.5*(Windspd[idx]+Windspd[idx+1])

    beta = beta_max * sall_temp_effect(T_used)  # pathogen growth rate

    dydt = np.zeros_like(y)
# =============================================================================
#     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     dydt[0] = %YOUR CODE GOES HERE for our Pb (Berries) function
#     dydt[1] = %YOUR CODE GOES HERE for our P function
#     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# =============================================================================

    # Disease variables
    dydt[2] = -beta * S * I + dydt[1] / A       # change in S
    dydt[3] = beta * S * I - mu_L_used * L + e  # change in L
    dydt[4] = mu_L_used * L - mu_I * I          # change in I
    dydt[5] = mu_I * I                          # change in R

# =============================================================================
#     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     %YOUR CODE GOES HERE for our E function
    if(I==0):  #spore production shouldn't start before infection (quirk of exponential curve fit)
        dydt[7] = 0
    else:
#       YOUR CODE GOES HERE for our F function
#     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# =============================================================================
    return dydt
