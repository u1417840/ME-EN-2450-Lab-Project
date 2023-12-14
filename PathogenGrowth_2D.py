# =============================================================================
# % function to compute the increase in icidence of a plant pathogen (or any
# % pathogen with stationary hosts).  This code is 2-D because it predicts 
# % the growth of the pathogen for a 2-D array of hosts each with a 
# % population that grows in size with time.  Additionally, it calculates 
# % airborne transmission of disease between hosts in the 2D array.  It is 
# % a modified version of the PathogenGrowth_0D function used in lab 09.
# %
# % The primary equations governing this are compound interest equations for
# % growth:
# %
# % dS/dt = -beta*S*I+dP/dt
# %
# % dL/dt = beta*S*I-mu_L*L+dE/dt
# %
# % dI/dt = mu_L*L*I-mu_I*I
# %
# % dP/dt = dB/dt + d(leaf)/dt 
# %
# % dR/dt = mu_I*I
# %
# % dE/dt = e
# %
# % dF/dt = Gamma*exp(alpha*I) - F*R_frac
# %
# % with:
# % S    = susceptible fraction of population (susceptible tissue)
# % beta = rate of infection growth for healthy population (fraction per day)
# % L    = fraction of tissue infected and latent (e.g., dormant before infection)
# % I    = fraction of tissue infected and producing inoculum
# % R    = fraction of tissue recovered (or removed) from population
# % P    = size of the total population (plant surface area)
# % B    = surface area of berries
# % E    = amount (fraction of population) of introduced disease from external sources
# % F    = size of the spreading population, e.g. sporulating for a fungus 
# % t    = time
# % mu_L = rate of decrease in latent population (fraction per day)
# % mu_I = rate of decrease in infectious population (fraction per day)
# % e    = rate of import from external sources
# %
# % inputs: vine (structure containing the initial size of susceptible population); beta; mu;
# % tspan (days to simulate array);
# % output: S,L,I,R,P,E,time (vector of simulation times), and B
# =============================================================================

import numpy as np
from SLIRPE import SLIRPE_model
from gaussian_plume_dep import gaussian_plume_dep
from latentperiod import latentperiod
# =============================================================================
# # IMPORT a function for time integration (can be Euler or RK4 or ...)
from euler import euler
# =============================================================================

def PathogenGrowth_2D(vine, beta_max, mu_L_target, mu_I, A, eta, kappa, xi, Gamma, alpha, T, U, V, tspan, NpX, NpY, Nsteps):

    # Set parameters in a list
    p = [beta_max, 1 / mu_I, T, tspan, A, np.sqrt(U**2 + V**2), np.degrees(np.arctan2(V, U)), eta, kappa, xi, Gamma, alpha]
   
    # Set parameters for scouting function
    ScoutSpeed = 0.125 # m/s
    DetectSize = 20 * ScoutSpeed # infection diameter (mm) that we are able to detect
    num_drones = 6
    scouts_per_day = 1
    drone_cost = 100 # cost of using drone for 1 hour to scout
    day_cost = 1000 # charge of not finding the infection each day
    day_cost_start = 10 # day that we start getting charged for not finding infection
    day_start_scout = int((day_cost_start + (mu_L_target / 2) - 0.5) * 24)
    drone_scout_cost = 0 # initial cost of using drones
    day_undetection_cost = 1500 # inital cost cause we start scouting on day 12.5
    scout_index = 1 # tracks if we found Infection, so we don't have to keep scouting the plants after we find it the first time
    
    def scout_fun(NpX, NpY, vine, A, DetectSize, num_drones, scouts_per_day, drone_scout_cost, day_undetection_cost, t):
        for i in range(NpX):
            for j in range(NpY):
                cnt = i + j * NpX
                if vine[cnt]["IsInfect"] == True:
                    I_dia = np.power( (vine[cnt]["I"][t] * (4 * A) / np.pi), 0.5)
                    if DetectSize <= (I_dia * 10):
                        drone_scout_cost = drone_scout_cost + (drone_cost * num_drones * scouts_per_day)
                        day_undetection_cost = day_undetection_cost + day_cost
                        return True, drone_scout_cost, day_undetection_cost
                    
        drone_scout_cost = drone_scout_cost + (drone_cost * num_drones * scouts_per_day)
        day_undetection_cost = day_undetection_cost + day_cost
        return False, drone_scout_cost, day_undetection_cost
    
    # Declare function handles
    odefun = lambda t, y, e, g: SLIRPE_model(t, y, e, g, p)

    # Loop over time steps
    for t in range(1,Nsteps):
        print(f"day={tspan[t]:.2f} infected plants={np.sum([v['IsInfect'] for v in vine])}")
        dt = tspan[t] - tspan[t - 1]  # timestep

        # Update list of infected vines
        active_vines = [idx for idx, v in enumerate(vine) if v['IsInfect']]
        
        # Initialize deposition flux for timestep
        dep_flux_sum = np.zeros(NpX * NpY)

        for idx in active_vines:
            # Calculate positions of other vines relative to each infected plant
            Xplume = [v["X"] - vine[idx]["X"] for v in vine]
            Yplume = [v["Y"] - vine[idx]["Y"] for v in vine]
          
            dep_flux = gaussian_plume_dep(Xplume, Yplume, p[5][t], p[6][t], vine[idx]["S"][t - 1], vine[idx]["F"][t - 1])  # /dt
            nnan_ind = ~np.isnan(dep_flux)
           
            dep_flux_sum[nnan_ind] += dep_flux[nnan_ind]
            
        # Loop over all vines
        for i in range(NpX):
            for j in range(NpY):
                cnt = i + j * NpX  # index counter for vectorized vine structure

                # Check if vines have just become latent
                if vine[cnt]["L"][t - 1] > 1e-8 and not vine[cnt]["LatentSwitch"]:
                    vine[cnt]["mu_L"] = latentperiod(t, dt, Nsteps, mu_L_target, np.zeros_like(T), T)
                    vine[cnt]["LatentSwitch"] = True

                # Set initial conditions for time integration
                y0 = [
                    vine[cnt]["B"][t - 1],  # amount of population, surface area, that is berries
                    vine[cnt]["P"][t - 1],  # total population, surface area, including berries and leaves
                    vine[cnt]["S"][t - 1],  # initial susceptible population fraction
                    vine[cnt]["L"][t - 1],  # initial latent population fraction
                    vine[cnt]["I"][t - 1],  # initial infectious population fraction
                    vine[cnt]["R"][t - 1],  # initial recovered population fraction
                    vine[cnt]["E"][t - 1],  # initial external population fraction
                    vine[cnt]["F"][t - 1],  # size of the spreading population
                ]

#             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#             % integrate using a time integration method from class
#             %%%% INSERT YOUR CODE HERE for time integration, note for
#             %%%% odefun to work as given above your call needs to look like:
#             %[y] = TimeInt(odefun,t,dt,y0,DepFlux_sum(cnt),vine(cnt).mu_L)
                # need to modify to match above^^
                y = euler(odefun, t, dt, y0, dep_flux_sum[cnt],vine[cnt]["mu_L"])
#             %NOTE: recognize that you are only integrating 1 time step!
#             %your routine can be more general than that but recognize that
#             %this point is in the middle of a time loop!
#             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                # Set outputs
                vine[cnt]["B"][t] = y[0]
                vine[cnt]["P"][t] = y[1]
                vine[cnt]["S"][t] = y[2]
                vine[cnt]["L"][t] = y[3]
                vine[cnt]["I"][t] = y[4]
                vine[cnt]["R"][t] = y[5]
                vine[cnt]["E"][t] = y[6]
                vine[cnt]["F"][t] = y[7]

                # Define a threshold for dispersal to start
                if vine[cnt]["I"][t] >= (np.pi * 0.25**2) / (4 *  A):
                    vine[cnt]["IsInfect"] = True

                # Turn off dispersal if we fall below the above size
                if vine[cnt]["I"][t] < (np.pi * 0.25**2) / (4 *  A) and vine[cnt]["IsInfect"]:
                    vine[cnt]["IsInfect"] = False


# =============================================================================
#     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     %%%        RECOMMENDED LOCATION FOR YOUR SCOUTING ROUTINE           %%%
#     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# =============================================================================
        if scout_index < 2: # checks to see if we still need to scout
                    for sctDay in range(int((len(tspan) - day_start_scout) / 24)): 
                        # scouts once every 24 hours
                        if t == day_start_scout + sctDay * 24: 
                            if scout_fun(NpX, NpY, vine, A, DetectSize, num_drones, scouts_per_day, drone_scout_cost, day_undetection_cost, t)[0] == True:
                                scout_index = scout_index + 1
                                day_discovered = tspan[t]
                                drone_scout_cost = scout_fun(NpX, NpY, vine, A, DetectSize, num_drones, scouts_per_day, drone_scout_cost, day_undetection_cost, t)[1]
                                day_undetection_cost = scout_fun(NpX, NpY, vine, A, DetectSize, num_drones, scouts_per_day, drone_scout_cost, day_undetection_cost, t)[2]
                                scouting_total_cost = drone_scout_cost + day_undetection_cost
                                print('\n  Day:', tspan[t], 'Infection Discovered')
                                print('    Drone Scouting Operation Cost: $', drone_cost * num_drones * scouts_per_day, '\n')
                            else:
                                drone_scout_cost = scout_fun(NpX, NpY, vine, A, DetectSize, num_drones, scouts_per_day, drone_scout_cost, day_undetection_cost, t)[1]
                                day_undetection_cost = scout_fun(NpX, NpY, vine, A, DetectSize, num_drones, scouts_per_day, drone_scout_cost, day_undetection_cost, t)[2]
                                print('\n  Day:', tspan[t], 'Infection Not Found')
                                print('    Drone Scouting Operation Cost: $', drone_cost * num_drones * scouts_per_day, '\n')
                             
                # Stops simulation once Infection is found so code doesn't run for infinite years
                if scout_index > 1:
                    break
        
    print('\n  Infection Found on Day', day_discovered)
    print('  Total Cost: $', scouting_total_cost)
    print('    Total Drone Scouting Cost: $', drone_scout_cost)
    print('    Total Delayed Detection Cost: $', day_undetection_cost)
    
    
