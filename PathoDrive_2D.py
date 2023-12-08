import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from PathogenGrowth_2D import PathogenGrowth_2D

# driver script for our 2 dimensional pathogen simulation

# load Environmental Forcing data (U, V, T, tspan)
mat_data = scipy.io.loadmat('EnvironmentalForcing.mat')
T = mat_data["T"][0]
tspan = mat_data["tspan"][0]
U = mat_data["U"][0]
V = mat_data["U"][0]

# declare global variables so the function has access without passing
#global NpX, NpY, Nsteps

# set simulation constants
# Global domain constants (apply everywhere)
NpX = 50  # number of plants in the X-direction
NpY = 50  # number of plants in the Y-direction
Nsteps = len(T)  # number of time steps for integration

# Global plant parameters
LpX = 1  # X-direction physical length of a plant in meters
LpY = 1  # Y-direction physical length of a plant in meters
# normalization factor for a plant ('final' susceptible plant surface area in cm^2)
A = 5000
# initial average size of an individual plant in the population
# (=model after 30 days with constant temp of 15C + 1 for initial berry size)
P_i_ave = 1.33 * 30 * (-0.35968 + 0.10789 * 15 - 0.00214 * 15 * 15) * 30 + 1
P_i_std = 0.2 * P_i_ave  # variance in the initial growth (fraction of the average)

# Global pathogen parameters
beta_max = 2  # max rate infection spread under ideal conditions (1/day)
mu_L_min = 6  # min length of latent period (min number of days latent)
mu_I = 10  # rate infection clears (number of days infectious)
eta = 1  # release fraction scale factor
kappa = 0.75  # release fraction scale factor
xi = -2.0  # release fraction offset
Gamma = 1e-2  # spore production from Calonnec et al 2009 approx scaled as surface area coverage
alpha = 0.314  # spore production 2nd factor

# Initialize individual Plants (vines)
# Here we will use a structure (vine) to store all the different variables
# to keep the association between variables and locations
vn = ['X', 'Y', 'IsInfect', 'LatentSwitch', 'P', 'B', 'S', 'L', 'I', 'R', 'E', 'F', 'mu_L']
# list of variable names to store in the structure
# initialize the storage structure
vine = [{'X': 0, 'Y': 0, 'IsInfect': False, 'LatentSwitch': False,
         'P': np.zeros(Nsteps), 'B': np.zeros(Nsteps), 'S': np.zeros(Nsteps),
         'L': np.zeros(Nsteps), 'I': np.zeros(Nsteps), 'R': np.zeros(Nsteps),
         'E': np.zeros(Nsteps), 'F': np.zeros(Nsteps), 'mu_L': np.zeros(Nsteps)}
        for _ in range(NpX * NpY)]
# set the position and initial variables.
# Here we are using a vector for later convenience when looping over only subsets of vines with active infections
for i in range(NpX):
    for j in range(NpY):
        cnt = i + j * NpX  # counter to vectorize vine
        vine[cnt]['X'] = i * LpX - 0.5  # x-position of the center of a vine in meters
        vine[cnt]['Y'] = j * LpY - 0.5  # y-position of the center of a vine in meters
        vine[cnt]['P'][0] = P_i_ave + P_i_std * np.random.randn()  # initial population size (random around the average)
        vine[cnt]['B'][0] = 1  # initial size of the berry population (assumed small (1cm^2))
        vine[cnt]['S'][0] = vine[cnt]['P'][0] / A  # initial size of the susceptible population (normalized)
        vine[cnt]['L'][0] = 0  # initial fraction of the population that is latent
        vine[cnt]['I'][0] = 0  # initial fraction of the population that is infectious
        vine[cnt]['R'][0] = vine[cnt]['I'][0] * mu_I  # initial fraction of the population that is recovered
        vine[cnt]['E'][0] = 0  # initial fraction of population from external sources
        vine[cnt]['F'][0] = 0  # initial amount of 'spores' for spreading

# Randomly select an initial plant for the infection to start at and give
# it the equivalent of 1 small (0.5cm diameter) spot on a leaf that is infected (latent)
RandV = np.random.randint(NpX * NpY)  # a random integer value in the range of the number of vines
vine[RandV]['L'][0] = (np.pi * 0.25**2)/(4*A)

# call the pathogen function
# tic
PathogenGrowth_2D(vine,beta_max,mu_L_min,mu_I,A,eta,kappa,xi,Gamma,alpha,T,U,V,tspan,NpX,NpY,Nsteps);
# toc

# Calculate stats for the entire domain at each time step
S_ave = np.mean([v['S'] for v in vine], axis=0)
L_ave = np.mean([v['L'] for v in vine], axis=0)
I_ave = np.mean([v['I'] for v in vine], axis=0)
R_ave = np.mean([v['R'] for v in vine], axis=0)
P_ave = np.mean([v['P'] for v in vine], axis=0)
E_ave = np.mean([v['E'] for v in vine], axis=0)
F_ave = np.mean([v['F'] for v in vine], axis=0)
B_ave = np.mean([v['B'] for v in vine], axis=0)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% plot results
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
#%INSERT YOUR CODE HERE to fill in the rest...
plt.figure()
plt.plot(tspan, P_ave / A, '-k', tspan, B_ave / A, '--k', tspan, S_ave, '-.m', tspan,
         L_ave, '--g', tspan, I_ave, ':b', tspan, R_ave, '-.r', tspan, E_ave, ':r',
         tspan, F_ave, '--y', linewidth=2)
plt.legend(['Total Population', 'Berry Population', 'Susceptible', 'Latent',
            'Infected', 'Removed', 'External', 'Spores'], loc='best')
plt.xlabel('time (days)')
plt.ylabel('Population (fraction of initial)')
plt.title('average epidemic')
plt.gca().set(xlim=[0, 61])
plt.grid(True)
plt.show()

# Plot the initially infected vine of the field
plt.figure()
plt.plot(tspan, vine[RandV]['P'] / A, '-k', tspan, vine[RandV]['B'] / A, '--k', tspan,
         vine[RandV]['S'], '-.m', tspan, vine[RandV]['L'], '--g', tspan, vine[RandV]['I'],
         ':b', tspan, vine[RandV]['R'], '-.r', tspan, vine[RandV]['E'], ':r',
         tspan, vine[RandV]['F'], '--y', linewidth=2)
plt.legend(['Total Population', 'Berry Population', 'Susceptible', 'Latent',
            'Infected', 'Removed', 'External', 'Spores'], loc='best')
plt.xlabel('time (days)')
plt.ylabel('Population (fraction of initial)')
plt.title(f'Initial infection at X={vine[RandV]["X"]}, Y={vine[RandV]["Y"]}')
plt.gca().set(xlim=[0, 61])
plt.grid(True)
plt.show()

#%INSERT YOUR CODE HERE to add plotting of other elements and optional
#%things like movies



