#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that compares the results from the spectral method (single domain), with
different spatial and temporal discretizations
@author: Antonio Preziosi-Ribero
CFD - Pontificia Universidad Javeriana
"""


# Importing the functions that calculate with different discretizations
from main_f import spectral1D


# Importing libraries to plot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style


# Declaring variables to be run and then plotted
dT = np.array([0.5, 0.2, 0.1, 0.01, 0.005, 0.001, 0.0005, 0.0001])
Nx = np.array([5, 8, 10, 20, 50, 75, 80])

# Matrix to store maximum errors of each run
MD = np.zeros((np.size(dT), np.size(Nx)))

# Looping to find maximum error of each run
for i in range(0, np.size(dT)):
    
    for j in range(0, np.size(Nx)):
        
        MD[i, j] = np.max(spectral1D(Nx[j], dT[i]))


# Plotting the error curves for each Sx
plt.figure(1, figsize=(11, 8.5))
style.use('ggplot')        

plt.subplot(2, 2, 1)
line1 = plt.semilogy(Nx, MD[0,:], 'b')
plt.xlim([np.min(Nx), np.max(Nx)])
plt.ylim([1e-5, 1e0])
plt.ylabel(r'Infinity error norm')
#plt.legend(loc=3)
plt.title('dT = ' + str(dT[0])) 
#
plt.subplot(2, 2, 2)
line1 = plt.semilogy(Nx, MD[1,:], 'b')
plt.xlim([np.min(Nx), np.max(Nx)])
plt.ylim([1e-5, 1e0])
#plt.legend(loc=3)
plt.title('dT = ' + str(dT[1]))
#
plt.subplot(2, 2, 3)
line1 = plt.semilogy(Nx, MD[2,:], 'b')
plt.xlim([np.min(Nx), np.max(Nx)])
plt.ylim([1e-5, 1e0])
plt.ylabel(r'Infinity error norm')
plt.xlabel('Number of nodes in discretization')
#plt.legend(loc=3)
plt.title('dT = ' + str(dT[2]))
#
plt.subplot(2, 2, 4)
line1 = plt.semilogy(Nx, MD[3,:], 'b')
plt.xlim([np.min(Nx), np.max(Nx)])
plt.ylim([1e-5, 1e0])
plt.xlabel('Number of nodes in discretization')
#plt.legend(loc=3)
plt.title('dT = ' + str(dT[3]))
plt.draw()

plt.suptitle('Discretization comparison')