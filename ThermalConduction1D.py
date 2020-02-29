### 1D Heat Conduction using FDM & Euler's

### See: https://www.youtube.com/watch?v=EGFfQt0 for reference.
### The objective of this script is to become familiar with using the Finite 
### difference method in conjuction with Euler's method to solve  the heat 
### equation. Using a 1-Dimensional model to simplify the problem and improve 
### understanding. 

### Also this will be used as a test for using dictionary's to support multiple
### materials I wish to implement in the 2D model.

import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns


global Nx, Nt, dx, dt


def solve_FDM(topology):
    '''
    Solves the heat equation using finite difference method.
    '''
    global Nx, Nt, dx, dt
    
    for tdx in range(Nt-1): #solves for time + 1, hence Nt-1 
        for idx in range(1,Nx-1): #boundary conditions are fixed, hence range(1,Nx-1)
            topology[tdx+1,idx]['T'] = topology[tdx,idx]['T'] + (topology[tdx,idx]['kT']*dt)/(dx**2)*(topology[tdx,idx+1]['T'] - 2*topology[tdx,idx]['T'] + topology[tdx,idx-1]['T'])
    return topology



def showTop(topology,prop):
    '''
    Returns a 2D array of a specified property across time series
    '''
    global Nx, Nt 
    
    visTop = np.zeros((Nt,Nx))

    for itt in range(Nt):
        for itX in range(Nx):         
            visTop[itt,itX] = topology[itt,itX][prop]
            
    return visTop


#-----------------------------Model Properties--------------------------------


total_length = 10                           # length [m]
dx = 1                                      # increment of length [m]
Nx =  int(total_length/dx)                  # number of X increments [-]

total_time = 10                             # total duratiom [s]
dt = 0.25                                   # increment of time [s]
Nt = int(total_time/dt)                     # number of t increments [-]


properties = {'T':0,'kT': 2} #Based on example used in video.

topology = [[properties.copy() for idx in range(Nx)] for tdx in range(Nt)]
topology = np.array(topology)

# Boundary Conditions
topology[:,:1] = {'T':200,'kT': 2} #! There is an error when using topology[:,:1]['T'] to assign B.C's
topology[:,9:] = {'T':150,'kT': 2}

#-------------------------------Solving Model--------------------------------
# Check stability 
stab = lambda dt,dx,kT: dt*2*kT/(dx**2)
stab_val = stab(dt,dx,topology[0,0]['kT'])
#!!! I will need to work out how to determine stability for a model based on 
#multiple materials, where there will be different values of kT

if  stab_val > 1:
    # Stability Condition Fail
    print('Stability condition not met')
    print(stab_val)
    
else:
    # Stability Condition Pass
    print('A solution will now be found')
    # Find the solution
    topology = solve_FDM(topology)
    
    
#---------------------------------Results-------------------------------------
    #Produces the results as a heatmap
    yaxis_labels = np.round(np.linspace(0,total_time,Nt+1),decimals = 2)  
    xaxis_labels = np.round(np.linspace(0,total_length,Nx),decimals = 0)
    

    plt.rcParams["figure.figsize"]=5,8
    ax = sns.heatmap(showTop(topology,'T'),cmap = "plasma",yticklabels= yaxis_labels,xticklabels= xaxis_labels,cbar_kws={'label':'Temperature ($^{o}$C)'},square = False)
    
    plt.xlabel('Length (m)')
    plt.ylabel('Time (s)')
    plt.title('Heat conduction along a metal pipe with time')
    plt.show()