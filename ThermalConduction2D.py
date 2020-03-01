# 2D Heat Conduction Model: WIP, writing base structure


#           ???:Future Changes           !!!:Assumption

import numpy as np
import copy
import seaborn as sns
from matplotlib import pyplot as plt


global Nx, Ny, Nt

#---------------------------------Functions-----------------------------------
def populateTopology(topology,material,T):
    '''
    Fills the model array (toplogy) with IVP conditions and material properties
    ???This function will eventually import data defining model properties and 
    IVP conditions. Or will create script to do-so. For now...
    '''
    global Nx, Ny, Nt 
  
    for itY in range(Ny):
        for itX in range(Nx):  
            
             topology[itY][itX] = material.copy() # assign material to cell
             
             #??? temporary temperature IVP replace with function paramater
             if itX == 0:
                  topology[itY][itX]['T'] = T
                  topology[itY][itX]['BC'] = True 
             if itX == Nx-1:
                  topology[itY][itX]['T'] = 40
                  topology[itY][itX]['BC'] = True  
       
    return topology


def showTop(topology,prop,t):
    '''
    Returns a 2D array of a specified property @ a specified time
    '''
    global Nx, Ny, Nt 
    
    visTop = np.zeros((Ny,Nx))
    for itY in range(Ny):
        for itX in range(Nx):
            
            visTop[itY,itX] = topology[t,itY,itX][prop]
            
    return visTop
    

#------------------------Material Physical Properties-------------------------


#Thermal Resitance [W/M^2.degC]: Coulson & Richardson's Chemical Engineering Volume 1
# Table 1.1: Thermal Conducitvity of Selected Materials. 
#!!! Assuming constant thermal resisitivy with temperature at this moment
    #??? possibility to replace with functions

iron = {'kT': 48} #??? add in the other material properties in future
steel = {'kT': 45}
bronze = {'kt':189}
water = {'kT':0.66}

#!!! To get the model up and running I am currently leaning on the parameters given
# in the matlab model given in the references.
# "Heat Diffusion in 2D Square Plate Using Finite Difference Method with Steady-State Solution" 
# by Amr Mousa Mohamed
silver = {'T':0,'kT':629,'cp':233,'rho':10490,'BC':False}




#-----------------------------Model Properties--------------------------------

total_length = 1            #length in [m]
total_width = 1           #width in [m]
dx = 0.025                   #increment size in X and Y directions [m]
dy= 0.025

Nx = int(total_length / dx)      #number of increments in X & Y directiosn [-]
Ny = int(total_width / dy)


total_time = 500        #final time condition [s]
dt = 0.6                  #time interval [s]
Nt = int(total_time / dt)    #number of increments in time

properties = {'T':0,'kT': 0,'cp':0,'rho':0,'BC':False} #dict of properties each cell will contain
                              #T: Temperature [degC]
                              #kT: Thermal Resistance [W/M.degC]
                              #cp: Heat Capacity (constant pressure) [J/kg.K]
                              #rho: Material Desnsity [kg/m^3]

                              #B.C: Boolean term to identify cells containing boundary condition



#------------------------------Creating Model--------------------------------
#creating array to store results                        
topology = [[properties.copy() for itX in range(Nx)] for itY in range(Ny)] #create 2D array

topology = populateTopology(topology,silver,75) #assign Boundary Conditions (B.Cs) and properties

# An array of objects when copied using .copy() will treat each 3rd axis array as one,
# therefore copy.deepcopy() is used to address this issue.
topology = np.array([copy.deepcopy(topology) for itt in range(Nt)]) #duplicate for time series


# Cell reference notation
# topology[3,1,1] -> {'T': 0, 'kT': 0}        indexing whole cell
# topology[3,1,1]['T'] -> 0                   indexing indv property/value





#-------------------------------Solving Model--------------------------------


# Equation Functions:
#     Function for collected term 'alpha'
alphaF = lambda top: top[0,0,0]['kT']/(top[0,0,0]['cp']*top[0,0,0]['rho'])
#     Function for FDM Euler approximation of the 2-D heat equation
heat_eqn = lambda top, alpha, dt, dx, dy, itt, itx, ity: top[itt,ity,itx]['T'] + dt*alpha*( ((top[itt,ity-1,itx]['T']-2*top[itt,ity,itx]['T']+top[itt,ity+1,itx]['T'])/dx**2) + ((top[itt,ity,itx-1]['T']-2*top[itt,ity,itx]['T']+top[itt,ity,itx+1]['T'])/dy**2) )



stab = lambda dt,dx,dy,alpha: dt*2*alpha*((dx**-2) + (dy**-2)) # stability condition value

alphaV = alphaF(topology) # assigning to variables for convinience, in the future it will simply be just the functions
stabV = stab(dt,dx,dy,alphaV)

stab_pass = False #Gatekeeper for stability check

# Checking model stability for solving
if stabV > 1:
    print('Stability condition is not met')
else:
    print('Stability condition is met, model will be solved')
    stab_pass = True
# FDM & Euler's Solution:
    for itt in range(Nt-1):
        for itY in range(1,Ny-1):
            for itX in range(1,Nx-1):
                if topology[itt,itX,itY]['BC'] == False:
                    topology[itt+1,itY,itX]['T'] = heat_eqn(topology,alphaV,dt,dx,dy,itt,itX,itY)
                    
                    
                    

#---------------------------------Results-------------------------------------

    #???Final version will show begin/end as sublots, possibility for animation?


if stab_pass == True:

    t_plot = total_time
        
    plt.rcParams["figure.figsize"]= 8,7 # this will need to be calculated
    #ax = sns.heatmap(showTop(topology,'T',0),cmap = "plasma",yticklabels= yaxis_labels,xticklabels= xaxis_labels,cbar_kws={'label':'Temperature ($^{o}$C)'},square = False)
    ax = sns.heatmap(showTop(topology,'T',t_plot-1),cmap = "plasma",cbar_kws={'label':'Temperature ($^{o}$C)'},square = False)
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.title('Heat conduction through metal sheet at t = %d seconds' %t_plot)
    
    plt.xticks(np.linspace(0,Nx,int(Nx/2)),np.round(np.linspace(0,total_width,int(Nx/2)),decimals = 2),rotation = 'vertical')
    plt.yticks(np.linspace(Ny,0,int(Ny/2)),np.round(np.linspace(0,total_length,int(Ny/2)),decimals = 2),rotation = 'horizontal')
    
    plt.show()
