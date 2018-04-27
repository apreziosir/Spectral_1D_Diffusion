#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 16:56:03 2018
Mapping natural coordinates to real coordinates
@author: apreziosir
"""

import numpy as np

def mapping(LGLP, X0, XF):
    
    xn = np.zeros(len(LGLP))
    dx = XF - X0
    
    for i in range(0, len(LGLP)):
        
        xn[i] = X0 + (LGLP[i] + 1) * (dx / 2)
        
    return xn