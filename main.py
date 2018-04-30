#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 09:36:54 2018

@author: apreziosir
"""
# Importing modules that are going to be used in the program

import numpy as np
import scipy.sparse as sp
import mapping as mp
import Lag_pder as lpd
import Analyt as AN
import matplotlib.pyplot as plt
from matplotlib import style

# ==============================================================================
# Physical variables
# ==============================================================================

# Initial x coordinate of the phenomena (m)
X0 = 0.

# Final x coordinate of the phenomena (m)
XF = 5.

# Initial time analyzed for the phenomena (s)
T0 = 1.

# Final time for the phenomena (s)
TF = 2.

# Diffusion coefficient in (m2/s)
Dx = 0.3

# Mass of tracer injected to the system
M = 10

# Injection point in the domain
xo = 2

# ==============================================================================
# Numerical parameters
# ==============================================================================

# Stability parameter Sx for diffusion when handling explicit methods
Sx = 0.1

# Crank-Nicholson ponderation factor theta for specific cases
# theta = 0.0 explicit
# theta = 1.0 implicit
theta = 1.

# Number of nodes in the domain
N = 30

# Calculating domain length
L = XF - X0

# ==============================================================================
# Start different calculations given the parameters from the model
# ==============================================================================

# Calculating the dt with the stability parameter
dT = 0.00005

# Defining number of timesteps
nT = int((TF - T0) / dT)

# Error saving
ert = np.zeros(int(nT))

# d(chi)/dx - Constant value that arose from mapping to natural coordinates
dchi_dx = 2 / (XF - X0)

# Calculating GLL points, as N - 2 points and the adding -1 and 1 as extreme 
# values (this is already in natural coordinates)
[LGLP, W] = np.polynomial.legendre.leggauss(N - 2)
LGLP = np.insert(LGLP, 0, -1)
LGLP = np.append(LGLP, 1)

# Calculating real coordinates of points from natural coordinates location
xn = mp.mapping(LGLP, X0, XF)

# P coefficient for calculations
P = dchi_dx ** 2 * dT * Dx

# ==============================================================================
# Building the differentiation matrix. I know I need the second derivative 
# matrix, but it is a fact that I can square the first derivative matrix to 
# obtain the second derivative matrix
# ==============================================================================

# Defining diffusion matrix for the problem (as sparse matrix)
# K = sp.lil_matrix((N, N))
K = np.zeros((N,N))

# Filling each row of the matrix qith the derivatives of the Lagrange 
# polynomials
for i in range(0, N):
    
    K[:, i] = lpd.Lag_pder(LGLP, i)
    
# Calculating second derivative as the matrix matrix multiplication
D2 = np.matmul(K, K)

# Freeing space - take out to prove it works
del(K)

# Final matrix for solving
K = P * D2 - np.identity(N)
K[0, :] = np.zeros(N)
K[N - 1, :] = np.zeros(N)
K[0, 0] = 1.
K[N -1, N - 1] = 1.

plt.spy(K)
plt.draw()
plt.pause(1.5)
plt.clf()
# ==============================================================================
# Initial condition and start up of other parameters
# ==============================================================================

# Generating initial condition
C = AN.difuana(M, XF - X0, Dx, xn, xo, T0)

C1 = np.zeros(N)
Cmax = np.max(C)

# Plotting initial condition
plt.ion()
plt.figure(1, figsize=(11, 8.5))
style.use('ggplot')

plt.subplot(1, 1, 1)
plt.plot(xn, C)
plt.title('Initial condition')
plt.xlabel(r'Distance $(m)$')
plt.ylabel(r'Concentration $ \frac{kg}{m} $')
plt.draw()
plt.pause(1.5)

# Entering time loop 

for t in range(1, nT):
    
    # Generating analytical solution
    Ca = AN.difuana(M, XF - X0, Dx, xn, xo, T0 + t * dT)
    
    # Setting up right hand side
    C[0] = -AN.difuana(M, L, Dx, X0, xo, T0 + t * dT)
    C[N - 1] = -AN.difuana(M, L, Dx, XF, xo, T0 + t * dT)
    
    # Solving system (matrix vector multiplication)
    C1 = -np.linalg.solve(K, C)
    
    # Estimating error
    err = np.abs(C1 - Ca)
    ert[t] = np.linalg.norm(err)
    
    # Plotting numerical solution and comparison with analytical
    plt.clf()
    
    plt.subplot(2, 2, 1)
    plt.plot(xn, C1, 'b')
    plt.xlim([X0, XF])
    plt.ylim([0, Cmax])
    plt.ylabel(r'Concentration $ \frac{kg}{m} $')
    plt.title('Numerical solution')
    
    plt.subplot(2, 2, 2)
    plt.plot(xn, Ca)
    plt.xlim([X0, XF])
    plt.ylim([0, Cmax])
    plt.title('Analytical solution')
    
    plt.subplot(2, 2, 3)
    plt.semilogy(xn, err)
    plt.xlim([X0, XF])
    plt.ylim([1e-8, 1e2])
    plt.ylabel('Absolute error')
    plt.title('Error')
    
    plt.subplot(2, 2, 4)
    plt.semilogy(np.linspace(T0, TF, nT), ert)
    plt.xlim([T0 - 0.2, TF + 0.2])
    plt.ylim([1e-8, 1e2])
    plt.title('Error evolution')
    
    plt.draw()
    titulo = 'Spectral elements method in single domain solution implicit'
    plt.suptitle(titulo)
    plt.pause(0.2)
    
    # Preparing for next timestep   
    C = C1