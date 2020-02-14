import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import numpy as np




# PUBLIC FUNCTIONS

def setup_figure(dimensions):
    fig = plt.figure(figsize=(15,15))
    if (dimensions == 3):
        ax1 = plt.subplot2grid((1,1), (0,0),projection='3d')
    elif  (dimensions == 2):
        ax1 = plt.subplot2grid((1,1), (0,0))
        
        
    return ax1


def PlotTheOrbitalDynamics(path_to_data,dimensions,ax):

        
    #Load data
    data = np.loadtxt(path_to_data)
    r = data[:,0]
    theta = data[:,1]
    phi = data[:,2]
    a = data[0,3]


    #Convert from Boyer-Lindquist to Cartesian Coordinates
    mm = np.sqrt(r**2 + a**2)
    x = mm * np.sin(theta)*np.cos(phi)
    y = mm * np.sin(theta)*np.sin(phi)
    z = r  * np.cos(theta)
    
    if (dimensions == 3):
        ax.plot(x,y,z)  
        ax.scatter(0,0,0, c='r') 
        
        #And make it pretty
        ax.tick_params(labelbottom=False,labeltop=False, labelleft=False)    

        
    elif  (dimensions == 2):
        ax.plot(x,y)
        ax.scatter(0,0, c='r') 

        
        
 
 






 
def PlotTheRay(path_to_data,dimensions,ax):

        
    #Load data
    data = np.loadtxt(path_to_data)
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]


    if (dimensions == 3):

        ax.plot(x,y,z) 
        ax.scatter(x[-1],y[-1],z[-1],c='C0')  
 
        
        #And make it pretty
        ax.tick_params(labelbottom=False,labeltop=False, labelleft=False)    

        
    elif  (dimensions == 2):

        ax.plot(x,y)
        ax.scatter(x[-1],y[-1], c='C0') 

        
        limit = 60
        ax.set_xlim(-limit, +limit)
        ax.set_ylim(-limit, +limit)



