#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 10:44:45 2018
Function that estimates the values of the derivative of the Lagrange polynomials
evaluated in each point of the grid
@author: Antonio Preziosi-Ribero
CFD - Pontificia Universidad Javeriana
"""

import numpy as np

def Lag_pder(x, I):
    
    # Polynomial degree for derivatives
    k = len(x)
    # Point that is analyzed (the row of the differentiation matrix)
    j = I
    # Vector that will store the values of the derivatives evaluated in the 
    # given point
    Lag_pder = np.zeros(len(x))
    xj = x[I]
    
    # Looping over vector elements - k elements
    for z in range(0, k):
        
        xa = x[z]
        add0 = 0
        temp1 = 0
        
        # What comes inside the sum - i for each polynomial
        for i in range(0, k):
            
            if i != j:
                
                add0 = 1 / (xj - x[i])
                prod0 = 1. 
                
                # What comes inside the productory - looping over m
                for m in range(0, k):
                    
                    if m != j and m != i:
                        
                        prod0 *= (xa - x[m]) / (xj - x[m])
                
                temp1 += add0 * prod0
                        
        Lag_pder[z] = temp1
        
    return Lag_pder