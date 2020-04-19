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

#method = 'Euler'
method = 'RK'

total_length = 1            #length in [m]
total_width = 1           #width in [m]
dx = 0.025                   #increment size in X and Y directions [m]
dy= 0.025

Nx = int(total_length / dx)      #number of increments in X & Y directiosn [-]
Ny = int(total_width / dy)


total_time = 10        #final time condition [s]
dt = 0.6                  #time interval [s]
Nt = int(total_time / dt)    #number of increments in time

properties = {'T':0,'kT': 0,'cp':0,'rho':0,'BC':False} #dict of properties each cell will contain
                              #T: Temperature [degC]
                              #kT: Thermal Resistance [W/M.degC]
                              #cp: Heat Capacity (constant pressure) [J/kg.K]
                              #rho: Material Desnsity [kg/m^3]

                              #B.C: Boolean term to identify cells containing boundary condition

ss_tol = 0.01 #tolerance for steady state calculation (difference between max/min values between iterations)
dyn_tol = 0.01 #tolerance for dynamic calculations (difference between max/min and s.s value)


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
heat_eqn = lambda top, alpha, dt, dx, dy, itt, itx, ity: top[itt,ity,itx]['T'] + dt*alpha*( ((top[itt,ity-1,itx]['T']-2*top[itt,ity,itx]['T']+top[itt,ity+1,itx]['T'])/dy**2) + ((top[itt,ity,itx-1]['T']-2*top[itt,ity,itx]['T']+top[itt,ity,itx+1]['T'])/dx**2) )


#     T''(x,y) for Runge-Kutta method
RK_eqn = lambda T, alpha, dx, dy, itx, ity: alpha*( ((T[itx-1,ity]-2*T[itx,ity]+T[itx+1,ity])/dy**2) + ((T[itx,ity-1]-2*T[itx,ity]+T[itx,ity+1])/dx**2))

#     FDM Heat equation re-arranged for steady state conditions where dT/dt = 0 
ss_heat_eqn = lambda top, dx, dy, itx, ity: (dx**2*(top[ity-1,itx]+top[ity+1,itx])+dy**2*(top[ity,itx-1]+top[ity,itx+1]))/(2*(dx**2 + dy**2))

#     Stability condition function
stab = lambda dt,dx,dy,alpha: dt*2*alpha*((dx**-2) + (dy**-2))




alphaV = alphaF(topology) # assigning to variables for convinience, in the future it will simply be just the functions
stabV = stab(dt,dx,dy,alphaV)  # stability constant of system

stab_pass = False #Gatekeeper for stability check

#Snapshot of initial temperatures for error checking later
IC_frame = showTop(topology,'T',0)
T_IC_max = np.max(np.max(IC_frame)) #greatest initial temperature
T_IC_min = np.min(np.min(IC_frame)) #smallest initial temperature

#-------------------------------~~Steady State~~------------------------------
topo_ss = showTop(topology,'T',0) # assign initial conditions from time = 0
# (steady-state values are independant of physical properties hence an array of temperatures is OK)

topo_ss_next = topo_ss.copy() # stores the next iteration of the calculation

max_ss_err = 1; min_ss_err = 1; # allocating error values
k = 1
while (max_ss_err >= ss_tol) or (min_ss_err >= ss_tol):
    #print(k)
    for ity in range(1,Ny-1):
        
        for itx in range(1,Nx-1):
            topo_ss_next[ity,itx] = ss_heat_eqn(topo_ss,dx,dy,itx,ity) # calculate T steady-state
            
    ss_diff = topo_ss_next-topo_ss
            
    max_ss_err = np.abs(np.max(np.max(ss_diff[np.nonzero(ss_diff)]))) # return the difference between the greatest values
    min_ss_err = np.abs(np.min(np.min(ss_diff[np.nonzero(ss_diff)]))) # return the difference between the greatest values
    
    topo_ss = topo_ss_next.copy()
    k += 1 # count to calculate time in a later update ( time to s.s = k * dt )
    


#-------------------------------~~Dynamic Section~~------------------------------
# Checking model stability for solving
    
if stabV > 1:
    print('Stability condition is not met')
else:
    print('Stability condition is met, model will be solved')
    stab_pass = True

# FDM & Euler's Solution:
    
    if method == 'Euler':
        
        euler_err_max = [] #creating error arrays
        euler_err_min = []
        
        for itt in range(Nt-1):
            
            if itt > 0:
                if (( 0 < euler_err_max[itt-1] <= dyn_tol) or (0 < euler_err_min[itt-1] <= dyn_tol)) :
                    print('Euler converged to tolerance within %d iterations' %(itt-1))
                    break     
                  
            for itY in range(1,Ny-1):
                
                for itX in range(1,Nx-1):
                    
                    if topology[itt,itX,itY]['BC'] == False:
                        
                        topology[itt+1,itY,itX]['T'] = heat_eqn(topology,alphaV,dt,dx,dy,itt,itX,itY)
            
            #Calculating and checking error
                        
            euler_frame = showTop(topology,'T',itt+1)
            euler_diff =  euler_frame - topo_ss
            
            euler_err_max.append(np.round(np.abs(np.max(np.max(euler_diff))),5))
            euler_err_min.append(np.round(np.abs(np.min(np.min(euler_diff))),5))
            
            if itt > 0:
                # Checking for incorrect convergance
                if ((euler_err_max[itt] == euler_err_max[itt-1] ) or (euler_err_min[itt] == euler_err_min[itt-1])) and (euler_err_max[itt] != 0):
                    print('Euler Error: Solution converging outside steady state')
                    stab_pass = False
                    break
            # Checking for an incorrect calculation if a temperautre is calculated outside the range of initial conditions    
            if (np.max(np.max(euler_frame)) > T_IC_max) or (np.min(np.min(euler_frame)) < T_IC_min):
                print('Euler Error: Solution not convering, temperature(s) outside initial conditions')
                stab_pass = False
                break
            
            
            
# 2nd Order Runge-Kutta                      
    elif method == 'RK':
        
        RK_err_max = [] # creating error arrays
        RK_err_min = []
        
        for itt in range(Nt-1):
            
            if itt > 0:
                if (( 0 < RK_err_max[itt-1] <= dyn_tol) or (0 < RK_err_min[itt-1] <= dyn_tol)) :
                    print('RK converged to tolerance within %d iterations' %(itt-1))
                    break     
            
            for itY in range(1,Ny-1):
                
                for itX in range(1,Nx-1):
                   
                    if topology[itt,itX,itY]['BC'] == False:                                         
                        
                        k1_T_vals = showTop(topology,'T',itt)
                     
                        k1 = RK_eqn(k1_T_vals, alphaV, dx, dy, itX, itY)
                                
                        k2_T_vals =  k1_T_vals + k1*dt
                                    
                        k2 = RK_eqn(k2_T_vals, alphaV, dx,dy ,itX ,itY) 
                        
                        topology[itt+1,itX,itY]['T'] = topology[itt,itX,itY]['T'] + (dt/2)*(k1 + k2)

            #Calculating and checking error
            
            RK_frame = showTop(topology,'T',itt+1)
            RK_diff =  RK_frame - topo_ss
            
            RK_err_max.append(np.round(np.abs(np.max(np.max(RK_diff))),5))
            RK_err_min.append(np.round(np.abs(np.min(np.min(RK_diff))),5))
            
            if itt > 0:
                # Checking for incorrect convergance
                if ((RK_err_max[itt] == RK_err_max[itt-1] ) or (RK_err_min[itt] == RK_err_min[itt-1])) and (RK_err_max[itt] != 0):
                    print('RK Error: Solution converging outside steady state')
                    stab_pass = False
                    break
            # Checking for an incorrect calculation if a temperautre is calculated outside the range of initial conditions
            if (np.max(np.max(RK_frame)) > T_IC_max) or (np.min(np.min(RK_frame)) < T_IC_min):
                print('RK Error: Solution not convering, temperature(s) outside initial conditions')
                stab_pass = False
                break

# Steady State Plot
                    
plt.rcParams["figure.figsize"]= 8,7 # this will need to be calculated
ax = sns.heatmap(topo_ss,cmap = "plasma",cbar_kws={'label':'Temperature ($^{o}$C)'},square = False)
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.title('Temperature Profile at Steady State')

plt.xticks(np.linspace(0,Nx,int(Nx/2)),np.round(np.linspace(0,total_width,int(Nx/2)),decimals = 2),rotation = 'vertical')
plt.yticks(np.linspace(Ny,0,int(Ny/2)),np.round(np.linspace(0,total_length,int(Ny/2)),decimals = 2),rotation = 'horizontal')

plt.show()


# Dynamic Plot

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
