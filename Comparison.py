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
dT = np.array([0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 
               0.00001])
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

plt.subplot(1, 2, 1)
for i in range(0, len(Nx)):
    line1 = plt.loglog(dT, MD[:, i], label= "N = " + str(Nx[i]))
#    plt.xlim([np.min(Nx), np.max(Nx)])
    plt.ylim([1e-6, 2e-1])
    plt.gca().invert_xaxis()
    plt.ylabel(r'Infinity error norm')
    plt.legend(loc=3)
    plt.title('Time refining error evolution') 

#
for i in range(0, len(dT)):
    plt.subplot(1, 2, 2)
    line1 = plt.semilogy(Nx, MD[i,:], label = 'dT = ' + str(dT[i]))
    plt.xlim([np.min(Nx), np.max(Nx)])
    plt.ylim([1e-6, 2e-1])
    plt.legend(loc=1)
    plt.title('Space refining error evolution')
#

plt.suptitle('Error analysis for spectral 1D diffusion')