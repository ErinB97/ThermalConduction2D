# 2D Heat Conduction Model: WIP, writing base structure


#           ???:Future Changes           !!!:Assumption

import numpy as np
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
            
             topology[itY][itX]['kT'] = material['kT'] #assign thermal resitance
           
             #temporary temperature IVP replace with function paramater
             if (itX < 15 or (itX > 85)) and ((itY < 15) or (itY > 135)):
                 topology[itY][itX]['T'] = T           
       
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
iron = {'kT': 48}
steel = {'kT': 45}
bronze = {'kt':189}
water = {'kT':0.66}


#-----------------------------Model Properties--------------------------------

total_length = 1            #length in [m]
total_width = 1.5           #width in [m]
dx = 0.01                   #increment size in X and Y directions [m]
dy= 0.01
dz = 0.1                    #depth of incrment size [m]

Nx = int(total_length / dx)      #number of increments in X & Y directiosn [-]
Ny = int(total_width / dy)


total_time = 500        #final time condition [s]
dt = 1                  #time interval [s]
Nt = int(total_time / dt)    #number of increments in time

properties = {'T':0,'kT': 0} #dict of properties each cell will contain
                              #T: Temperature [degC]
                              #kT: Thermal Resistance [W/M^2.degC]




#------------------------------Creating Model--------------------------------
#creating array to store results                        
topology = [[properties.copy() for itX in range(Nx)] for itY in range(Ny)] #create 2D array

topology = populateTopology(topology,steel,75) #assign IVP and properties

topology = [topology for itT in range(Nt)] #duplicate for time series
topology = np.array(topology) #convert to numpy array

# topology[3,1,1] -> {'T': 0, 'kT': 0}        indexing whole cell
# topology[3,1,1]['T'] -> 0                   indexing indv property/value




#-------------------------------Solving Model--------------------------------

    #???To be completed

#---------------------------------Results-------------------------------------

    #???Final version will show begin/end as sublots, possibility for animation?

fig, ax = plt.subplots()
im = ax.imshow(showTop(topology,'T'))
plt.show()