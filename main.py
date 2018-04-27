#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 09:36:54 2018

@author: apreziosir
"""
# Importing modules that are going to be used in the program

import numpy as np		# Python numerical library
import mapping as mp

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
TF = 5.

# Diffusion coefficient in (m2/s)
Dx = 0.3

# ==============================================================================
# Numerical parameters
# ==============================================================================

# Stability parameter Sx for diffusion when handling explicit methods
Sx = 0.

# Crank-Nicholson ponderation factor theta for specific cases
# theta = 0.0 explicit
# theta = 1.0 implicit
theta = 1.

# Number of nodes in the domain
N = 10

# ==============================================================================
# Start different calculations given the parameters from the model
# ==============================================================================

# Calculating the dt with the stability parameter
dt = Sx * (XF - X0) / Dx

# d(chi)/dx - Constant value that arose from mapping to natural coordinates
dchi_dx = 2 / (XF - X0)

# Calculating GLL points, as N - 2 points and the adding -1 and 1 as extreme 
# values (this is already in natural coordinates)
[LGLP, W] = np.polynomial.legendre.leggauss(N - 2)
LGLP = np.insert(LGLP, 0, -1)
LGLP = np.append(LGLP, 1)

# Calculating real coordinates of points from natural coordinates location
xn = mp.mapping(LGLP, X0, XF)

# ==============================================================================
# Building the differentiation matrix. I know I need the second derivative 
# matrix, but it is a fact that I can square the first derivative matrix to 
# obtain the second derivative matrix
# ==============================================================================

