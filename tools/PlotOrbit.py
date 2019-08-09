from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import os





d = 2


#Set up plotting environment
if (d == 3):
    fig = plt.figure(figsize=(10,10))
    ax1 = plt.subplot2grid((1,1), (0,0),projection='3d')
elif  (d == 2):
    fig = plt.figure(figsize=(20,10))
    ax1 = plt.subplot2grid((1,2), (0,0))
    ax2 = plt.subplot2grid((1,2), (0,1))



#Load data

path = os.environ['IliadDir']
data = np.loadtxt(path + 'MPDFormatData.txt')

r = data[:,0]
theta = data[:,1]
phi = data[:,2]
a = data[0,3]


mm = np.sqrt(r**2 + a**2)
x = mm * np.sin(theta)*np.cos(phi)
y = mm * np.sin(theta)*np.sin(phi)
z = r  * np.cos(theta)






#Plot it


if (d == 3):
    ax1.plot(x,y,z)  
    limit = max(max(x),max(y),max(z))
    ax1.set_xlim(-limit,+limit)
    ax1.set_ylim(-limit,+limit)
    ax1.set_zlim(-limit,+limit)

if (d == 2):
    ax1.plot(x,y)
    ax2.plot(x,z)






plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20

plt.show()
