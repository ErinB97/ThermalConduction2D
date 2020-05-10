# 2D Heat Conduction Model: 


import numpy as np
import copy
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd

# Importing non-mathematical functions from model_functions.py
from model_functions import *
from tkinter import simpledialog , filedialog



#------------------------Material Physical Properties-------------------------

#!!! To get the model up and running I am currently using the parameters given
# in the matlab model given in the references.
# "Heat Diffusion in 2D Square Plate Using Finite Difference Method with Steady-State Solution" 
# by Amr Mousa Mohamed

# An example for how Silver is stored in variable space
# silver = {'T':0,'kT':629,'cp':233,'rho':10490, 'alpha':2.57E-4,'BC':False}
                              # T: Temperature [degC]
                              # kT: Thermal Resistance [W/M.degC]
                              # cp: Heat Capacity (constant pressure) [J/kg.K]
                              # rho: Material Desnsity [kg/m^3]
                              # alpha: combined term for kT, cp, and density (kT/cp*density)
                              # B.C: Boolean term to identify cells containing boundary condition



#---------------------------Mathematical Functions-----------------------------


# Equation Functions:
#     Function for collected term 'alpha'
alphaF = lambda kT,cp,rho: kT/(cp*rho)

#     Function for FDM Euler approximation of the 2-D heat equation
euler_heat_eqn = lambda top, alpha,dt, dx, dy, itt, itx, ity: top[itt,ity,itx]['T'] + dt*alpha*( ((top[itt,ity-1,itx]['T']-2*top[itt,ity,itx]['T']+top[itt,ity+1,itx]['T'])/dy**2) + ((top[itt,ity,itx-1]['T']-2*top[itt,ity,itx]['T']+top[itt,ity,itx+1]['T'])/dx**2) )


#     Heat Equation function for Runge-Kutta method
RK_eqn = lambda T, alpha, dx, dy, itx, ity: alpha*( ((T[ity-1,itx]-2*T[ity,itx]+T[ity+1,itx])/dy**2) + ((T[ity,itx-1]-2*T[ity,itx]+T[ity,itx+1])/dx**2))

#     FDM Heat equation re-arranged for steady state conditions where dT/dt = 0 
ss_heat_eqn = lambda top, dx, dy, itx, ity: (dx**2*(top[ity-1,itx]+top[ity+1,itx])+dy**2*(top[ity,itx-1]+top[ity,itx+1]))/(2*(dx**2 + dy**2))

#     Stability condition function
stab = lambda dt,dx,dy,alpha: dt*2*alpha*((dx**-2) + (dy**-2))




#-----------------------Loading Model & Material Data--------------------------

# "L:1,W:1,dX:0.025,dY:0.025,dT:0.6,TT:10,ss_tol:0.01,dyn_tol:0.01"

# Loading material data:

material_path = filedialog.askopenfilenames(title = 'Select material data file path')[0]
material_file = pd.read_csv(material_path, index_col = 0)

# Loading model data:

model_path = filedialog.askopenfilenames(title = 'Select model file path')[0]
model_file = pd.read_csv(model_path, header = 0)


    # Obtain model properties from file header
total_length, total_width, Nx ,Ny, dx, dy, total_time, dt, Nt,  ss_tol, dyn_tol = GetModelParamaters(model_file.columns[0],model_file.shape) 

N = (Nt,Ny,Nx) # storing N values in tuple for passing to functions

# Creating and populating topology, and obtaining list of materials which are present.
topology, materials_used = populateTopology(model_file,material_file,alphaF,N)

# Cell reference notation
# topology[3,1,1] -> {'T': 0, 'kT': 0}        indexing whole cell
# topology[3,1,1]['T'] -> 0                   indexing indv property/value





#-------------------------------Model Set-Up----------------------------------


# Prompt user for which caculation method to use
method =  simpledialog.askstring('Calculation method', "Enter either 'euler' or 'RK' :")

# Checking the stability constant for all materials in the model

stab_pass = True #Gatekeeper for stability check

for mat in materials_used:
    # Calculate alpha for current material
    alphaV = alphaF(material_file.loc[mat,'kT'], material_file.loc[mat,'cp'], material_file.loc[mat,'density']) #value of alpha for the material
    stabV = stab(dt,dx,dy,alphaV)  # stability constant of material

    # check stability constant
    if stabV >  1:
        print('Stability condition for ',mat,' is not met. Failed at ',np.round(stabV,4),'. Increase dX and/or dY interval size. Or decrease dt.')
        stab_pass = False 
        

# The stability constant value check could be done within populateToplogy() because
# alpha is calculated for each cell. However, I wanted to do it at a higher level
# because this is a model calculation. 
#     
# Calculating alpha is an exception to this as each cell needs to have alpha 
# calculated from its properties, so is done within populateTopology()



#Storing initial temperatures for error checking later
IC_frame = showTop(topology,'T',0,N)
T_IC_max = np.max(np.max(IC_frame)) #greatest initial temperature
T_IC_min = np.min(np.min(IC_frame)) #smallest initial temperature




#-------------------------------Model Calculations-----------------------------

#-------------------------------~~Steady State~~------------------------------
topo_ss = showTop(topology,'T',0,N) # assign initial conditions from time = 0
# (steady-state values are independant of physical properties hence an array of temperatures is OK)

topo_ss_next = topo_ss.copy() # stores the next iteration of the calculation

max_ss_err = 1; min_ss_err = 1; # allocating error values

while (max_ss_err >= ss_tol) or (min_ss_err >= ss_tol):

    for ity in range(1,Ny-1):
        
        for itx in range(1,Nx-1):
            
            if topology[0,ity,itx]['BC'] == False:
                
                topo_ss_next[ity,itx] = ss_heat_eqn(topo_ss,dx,dy,itx,ity) # calculate T steady-state
            
    ss_diff = topo_ss_next-topo_ss
            
    max_ss_err = np.abs(np.max(np.max(ss_diff[np.nonzero(ss_diff)]))) # return the difference between the greatest values
    min_ss_err = np.abs(np.min(np.min(ss_diff[np.nonzero(ss_diff)]))) # return the difference between the greatest values
    
    topo_ss = topo_ss_next.copy() # making the next frame become the current frame, for the next calclulation step.
    



#-------------------------------~~Dynamic State~~------------------------------
    
# Checking model stability for solving
if stab_pass == True:
    print('Stability condition is met, model will be solved using %s method' %method)
    solve_time = total_time # assuming steady state is not reached, if SS is reached, will return time to SS.
    
    showConvergeWarning = False
    
# FDM & Euler's Solution:
    
   
    if method == 'euler':
        
        euler_err_max = [] #creating list to store error
        
        for itt in range(Nt-1):
            
            if itt > 0:
                # checking if calculation has converged to steady state
                if ( 0 < euler_err_max[itt-1] <= dyn_tol) :
                    # calculating time to reach steady state
                    
                    solve_time = int((itt-1)*dt)
                    print('Euler converged to S.S tolerance within %d iterations with a time of %d seconds' %(itt-1,solve_time))
                    print('Euler greatest difference to S.S values is: %f' %euler_err_max[itt-1])
                    break     
                    
            
            for itY in range(1,Ny-1):
                
                for itX in range(1,Nx-1):
                    
                    if topology[itt,itY,itX]['BC'] == False:
                        # calculate the temperature for cell @ time + dt
                        topology[itt+1,itY,itX]['T'] = euler_heat_eqn(topology,topology[itt,itY,itX]['alpha'],dt,dx,dy,itt,itX,itY)
            

            
            #Calculating and checking error
                        
            euler_frame = showTop(topology,'T',itt+1,N)
            euler_diff =  np.abs((euler_frame - topo_ss))
            
            euler_err_max.append(np.round(np.max(np.max(euler_diff[np.nonzero(euler_diff)])),3))
           
            
            if itt > 0:
                # Checking for incorrect convergance where error converges at non-zero value
                if (euler_err_max[itt] == euler_err_max[itt-1] )  and (euler_err_max[itt] != 0):                    
                    
                    showConvergeWarning = True
                    break
                    
                    
                    
            # Checking for an incorrect calculation if a temperautre is calculated outside the range of initial conditions    
            if (np.max(np.max(euler_frame)) > T_IC_max) or (np.min(np.min(euler_frame)) < T_IC_min):
                print('Euler Error: Solution not convering, temperature(s) outside initial conditions')
                stab_pass = False
                break
        # Give warning on possible bad convergance.
            
        if showConvergeWarning == True:
            solve_time = int((itt-1)*dt)
            print('Euler converged to 3 d.p at %.3f within %d iterations with a time of %d seconds' %(euler_err_max[-1],itt-1,solve_time))
            
        
        
# 2nd Order Runge-Kutta                      
    elif method == 'RK':
        
        RK_err_max = [] # creating list to store error
        
        for itt in range(Nt-1):
            
            if itt > 0:
                if ( 0 < RK_err_max[itt-1] <= dyn_tol)  :
                    solve_time = int((itt-1)*dt)
                    print('RK converged to S.S tolerance tolerance within %d iterations with a time of %d seconds' %(itt-1,solve_time))
                    print('RK greatest difference to S.S values is: %f' %RK_err_max[itt-1])
                    break     
            
            for itY in range(1,Ny-1):
                
                for itX in range(1,Nx-1):
                   
                    if topology[itt,itY,itX]['BC'] == False: 
                        
                        # calculate the temperature for cell @ time + dt                                        
                        
                        k1_T_vals = showTop(topology,'T',itt,N)
                     
                        k1 = RK_eqn(k1_T_vals, topology[itt,itY,itX]['alpha'], dx, dy, itX, itY)
                                
                        k2_T_vals =  k1_T_vals + k1*dt
                                    
                        k2 = RK_eqn(k2_T_vals, topology[itt,itY,itX]['alpha'], dx,dy ,itX ,itY) 
                        
                        topology[itt+1,itY,itX]['T'] = topology[itt,itY,itX]['T'] + (dt/2)*(k1 + k2)

            #Calculating and checking error
            
            RK_frame = showTop(topology,'T',itt+1,N)
            RK_diff = np.abs(RK_frame - topo_ss)
            
            RK_err_max.append(np.round(np.abs(np.max(np.max(RK_diff[np.nonzero(RK_diff)]))),3))
            
            
            if itt > 0:
                # Checking for incorrect convergance, where error converges at non-zero value
                if (RK_err_max[itt] == RK_err_max[itt-1] ) and (RK_err_max[itt] != 0):
                    showConvergeWarning = True
                    
            # Checking for an incorrect calculation if a temperautre is calculated outside the range of initial conditions
            if (np.max(np.max(RK_frame)) > T_IC_max) or (np.min(np.min(RK_frame)) < T_IC_min):
                print('RK Error: Solution not convering, temperature(s) outside initial conditions')
                stab_pass = False
                break
         
        # Give warning on possible bad convergance.
        if showConvergeWarning == True:
            solve_time = int((itt-1)*dt)
            print('RK converged to 3 d.p at %.3f within %d iterations with a time of %d seconds' %(RK_err_max[-1],itt-1,solve_time))

#----------------------------------Plotting Results---------------------------

dim = (total_width , total_length ) # size is needed for axis labels


# Dynamic Plot:

if stab_pass == True:
    
    plot_time = int(solve_time/dt)
    
    #Plot boundary condition
    ax = plotTopology(dim,N,topology,'BC','dyn', solve_time ,plot_time)
    plt.show()
    
    # Plot density
    ax = plotTopology(dim,N,topology,'rho','dyn', solve_time ,plot_time)
    plt.show()
    
    # Plot specific heat
    ax = plotTopology(dim,N,topology,'cp','dyn', solve_time ,plot_time)
    plt.show()
    
    # Plot thermal conductivity
    ax = plotTopology(dim,N,topology,'kT','dyn',solve_time ,plot_time)
    plt.show()

    
    # Plot temperature
    ax = plotTopology(dim,N,topology,'T','dyn',solve_time ,plot_time)
    plt.show()
    
    ax = plotTopology(dim,N,topology,'T','dyn',25 ,int(25/dt))
    plt.show()
    
    ax = plotTopology(dim,N,topology,'T','dyn')
    plt.show()
    
    # Plot error
    ax = plt.plot(euler_err_max)
    plt.xlabel('Iterations')
    plt.ylabel('Calculation Error')
    plt.title('Convergenace of Error ($^{o}$C)')
    plt.show()
    
# Steady State Plot:
               
ax = plotTopology(dim,N,topo_ss,'T','SS')
plt.show()
