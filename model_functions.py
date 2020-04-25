# Thermal Conduction Function File:

# This file will contain all the non-mathematic functions used by th Thermal 
# Conduction 2D model. 
#
# Mathematical functions and routines will be written in the main
# 'ThermalConduction2D file.

#import pandas as pd
#import seaborn as sns
#import numpy as np
#import copy
#from matplotlib import pyplot as plt

#global Nx, Ny, Nt, materials


def GetModelParamaters(input_str):
    '''
    This function will return model paramaters from a string where parameters 
    are deliminated by ',' and the key & value are separated by ':'
    
    Paramaters are returned as a dictionary
    '''
    
    input_str = input_str.split(',')
    
    parameters = {}
    
    for par in input_str:
        split_par = par.split(':')
        parameters.update({split_par[0]:float(split_par[1])})

    #For personal reference, notation being used in CSV file
    # L:1,W:1,dX:0.025,dY:0.025,dT:0.6,TT:10,ss_tol:0.01,dyn_tol:0.01

    total_length = parameters['L']            #length in [m]
    total_width =  parameters['W']            #width in [m]
    dx =  parameters['dX']                   #increment size in X and Y directions [m]
    dy=  parameters['dY'] 

    total_time = parameters['TT']        #final time condition [s]
    dt = parameters['dT']                 #time interval [s]
    
    ss_tol = parameters['ss_tol']
    dyn_tol = parameters['dyn_tol']
    
    return total_length, total_width, dx, dy, total_time, dt, ss_tol, dyn_tol




def populateTopology(file):
    '''
    Fills the model array (toplogy) with IVP conditions and material properties
    ???This function will eventually import data defining model properties and 
    IVP conditions. Or will create script to do-so. For now...
    '''
    global Nx, Ny, Nt 
    
    # Create the 2D array
    
    properties_temp = {'T':0,'kT': 0,'cp':0,'rho':0,'BC':False} # template dict of properties each cell will contain
    
    topology = [[properties_temp.copy() for itX in range(Nx)] for itY in range(Ny)] #create 2D array, list of lists of dicts
    
    
  
    dimensions = file.shape
    
    if (dimensions[0] == Ny) and (dimensions[1] == Nx): 
        # check the dimensions of the model in header align with actual dimensions given

        for itY in range(Ny):
            for itX in range(Nx):
                
                cell_str = file.iloc[itY,itX]
                cell_str = cell_str.split(',')
                
                # Fill the cell with values given in CSV cell
                topology[itY][itX] = populateCell(topology[itY][itX], cell_str)  
         
            
            
        topology = np.array([copy.deepcopy(topology) for itt in range(Nt)]) #duplicate for time series
        
        
    else:
        # !!! if the dimensions do not match, an error will be returned
        topology = False
       
    return topology




def populateCell(cell,cell_str):
    '''
    Fills cell given by populateTopology
    cell : dict
        cell of properties from topology
    cell_str : list of str
        list of str representing contents of cell from CSV file
    material : str
        reference of which material is in cell

    ''' 
    global materials
        
    # This has been separated out from populateTopology due to difference in scope
    
    properties = {}
    
    for prop in cell_str:
        prop_split = prop.split(':')
        properties.update({prop_split[0]:prop_split[1]})
        
    # T:0,M:silver,BC:False reference for CSV format
    # {'T':0,'kT':629,'cp':233,'rho':10490,'BC':False} reference for topology format
    
    #Assign cell properties   
      
    cell['T'] = properties['T']
    
    if properties['BC'] == 'True': # assign appropriate boolean value
        cell['BC'] == True
    else:
        cell['BC'] == False
    
    cell['kT'] = materials.loc[properties['M'],'kT']
    cell['cp'] = materials.loc[properties['M'],'cp']
    cell['rho'] = materials.loc[properties['M'],'density']


    return cell




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

def add(a,b):
    return a + b

