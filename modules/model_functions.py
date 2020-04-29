# Thermal Conduction Function File:

# This file will contain all the non-mathematic functions used by th Thermal 
# Conduction 2D model. 
#
# Mathematical functions and routines will be written in the main
# 'ThermalConduction2D file.

#import pandas as pd
import seaborn as sns
import numpy as np
import copy
from matplotlib import pyplot as plt

#global Nx, Ny, Nt, materials


def GetModelParamaters(input_str,shape):
    '''
    This function will return model paramaters from a string where parameters 
    are deliminated by ',' and the key & value are separated by ':'    
    '''
    
    input_str = input_str.split(',')
    
    parameters = {}
    
    for par in input_str:
        split_par = par.split(':')
        parameters.update({split_par[0]:float(split_par[1])})

    #For reference, notation being used in CSV file
    # L:1,W:1,dX:0.025,dY:0.025,dT:0.6,TT:10,ss_tol:0.01,dyn_tol:0.01

    total_length = parameters['L']       # length in [m]
    total_width =  parameters['W']       # width in [m]
    
    Ny, Nx = shape                       # number of Y and X intervals, given by number of cells
    
    dy = total_length / Ny               # interval size of Y and X intervals
    dx = total_width / Nx


    total_time = parameters['TT']        # final time condition [s]
    dt = parameters['dT']                # time interval [s]
    Nt = int(total_time / dt)            # number of time intervals
    
    ss_tol = parameters['ss_tol']        # tolerance for steady state calculation (difference between max/min values between iterations)
    dyn_tol = parameters['dyn_tol']      # tolerance for dynamic calculations (difference between max/min and s.s value)
    
 
    
    return total_length, total_width, Nx ,Ny, dx, dy, total_time, dt, Nt,  ss_tol, dyn_tol




def populateTopology(file,materials,alphaF,N):
    '''
    Fills the model array (toplogy) with IVP conditions and material properties,
    then duplicates across the time series.
    '''
    
    #alphaF = lambda kT,cp,rho: kT/(cp*rho)
    
    Nt, Ny, Nx = N
    
    # In order to check the stability constant of the whole model, need to know
    # which materials are present as the stab. constant is calculated from
    # physical properties
    materials_used = []
    
    # Create the 2D array
    
    properties_temp = {'T':0,'kT': 0,'cp':0,'rho':0,'alpha':0,'BC':False} # template dict of properties each cell will contain
    
    topology = [[properties_temp.copy() for itX in range(Nx)] for itY in range(Ny)] #create 2D array, list of lists of dicts    

    for itY in range(Ny):
        for itX in range(Nx):
            
            cell_str = file.iloc[itY,itX] # read string of cell properties from CSV cell
            cell_str = cell_str.split(',')
            
            # Fill the cell with values given in CSV cell
            topology[itY][itX], cell_mat = populateCell(topology[itY][itX], cell_str,materials,alphaF)  
          
            
            if materials_used.count(cell_mat) == 0:
                materials_used.append(cell_mat)
            
    topology = np.array([copy.deepcopy(topology) for itt in range(Nt)]) #duplicate for time series
        
        

       
    return topology, materials_used




def populateCell(cell,cell_str,materials,alphaF):
    '''
    Fills a cell's values as define by cell_str, called by populateTopology()
    ''' 
    
    properties = {}
    
    for prop in cell_str:
        prop_split = prop.split(':')
        properties.update({prop_split[0]:prop_split[1]})
        
    # T:0,M:silver,BC:False reference for CSV format
    # {'T':0,'kT':629,'cp':233,'rho':10490,'BC':False} reference for topology format
    
    #Assign cell properties:  
      
    cell['T'] = float(properties['T'])
    
    if properties['BC'] == 'True': # assign appropriate boolean value
        cell['BC'] = True
    else:
        cell['BC'] = False
    
    # Getting cell properties from materials table
    cell['kT'] = materials.loc[properties['M'],'kT']
    cell['cp'] = materials.loc[properties['M'],'cp']
    cell['rho'] = materials.loc[properties['M'],'density']

    cell['alpha'] = alphaF(cell['kT'],cell['cp'],cell['rho']) # in future round this 4 to sig.fig
    
    cell_mat = properties['M'] # used to check stability constant for all materials

    return cell,cell_mat



def showTop(topology,prop,itt,N):
    
    '''
    Returns a 2D array of a specified property @ a specified time
    '''
    #N = (Nt,Ny,Nx)
    Nt, Ny, Nx = N
    
    visTop = np.zeros((Ny,Nx))
    
    for itY in range(Ny):
        for itX in range(Nx):
            
            visTop[itY,itX] = topology[itt,itY,itX][prop]
      
    
    return visTop


def plotTopology(dim,N,topology,prop,mode,t_plot = 0 ,itt = 0):
    '''
    Returns pyplot axes object for the given property within a topology array.
    '''

    Nt,Ny,Nx = N
    total_width , total_length = dim
    
    # Extracts requested property from given topology
    if mode == 'dyn':
        if itt > 0:
            topology = showTop(topology, prop, itt-1, N)
        else:
            topology = showTop(topology, prop, itt, N)            
        
    # Strings for properties and their units for figure text
        # Displaying boundary condition to be added in future.
    fig_prop = {'T':'Temperature','cp':'Specific Heat','rho':'Density','kT':'Heat Conductivity','BC':'Boundary Condition'}
    fig_unit = {'T':'($^{o}$C)','cp':'(J / kg*K)','rho':'(kg / m$^{3}$)','kT':'(W/ m*K)','BC':'(-)'}
    
    # Create a title relevant to the type of calculation being plotted.
    if mode == 'dyn':
        title_str = ('Model at t = %d seconds showing %s profile' %(t_plot,fig_prop[prop]))
    elif mode == 'SS':
        title_str = ('%s profile at Steady State' %fig_prop[prop])
    
    
    # Creating the plot
    plt.rcParams["figure.figsize"]= 8,7
    
    ax = sns.heatmap(topology,cmap = "plasma",cbar_kws={'label':fig_prop[prop]+ ' ' + fig_unit[prop]},square = False)      
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')  

    # Creating tick marks that correspond to the dimensions of the model
    plt.xticks(np.linspace(0,Nx,int(Nx/2)),np.round(np.linspace(0,total_width,int(Nx/2)),decimals = 2),rotation = 'vertical')
    plt.yticks(np.linspace(Ny,0,int(Ny/2)),np.round(np.linspace(0,total_length,int(Ny/2)),decimals = 2),rotation = 'horizontal')
    
    plt.title(title_str)
    
    return ax