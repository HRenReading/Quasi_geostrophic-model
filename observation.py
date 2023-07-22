# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 20:06:05 2022

@author: 44754
"""


import os
import numpy as np
import numpy.random as rnd
import random
import scipy
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.io import FortranFile
from scipy.linalg import sqrtm
from scipy.linalg import circulant
import math


def vec2field(x,nxx,nyy):
    "function that transfers the vector to field"
    s1 = np.reshape(x[0:nxx*nyy],(nyy,nxx))
    s2 = np.reshape(x[nxx*nyy:2*nxx*nyy],(nyy,nxx))
    
    return s1,s2 


def measure(r,x,rng,Nx,L,nxx,nyy):
    "function maps the observations into observation space"
    Lnx = int(nxx-4)
    Lny = int(nyy)
    meas = np.zeros(2*Lny*Lnx)
    yrandom = r*rng.normal(0,1,2*Lny*Lnx)
    s1,s2 = vec2field(x,nxx,nyy)
    for j in range(Lny):
        for i in range(Lnx):
            meas[j*Lnx+i] = s1[j,i+2]
            meas[j*Lnx+i+Lnx*Lny] = s2[j,i+2]
            del i
        del j
    meas = meas + yrandom
    
    return meas
"""

def measure(r,x,rng,Nx,L,nxx,nyy):
    "function maps the observations into observation space"
    theta = 5
    meas = np.zeros(2*L)
    yrandom = r*rng.normal(0,1,2*L)
    x1 = x[0:nxx*nyy]
    x2 = x[nxx*nyy:]
    for l in range(L):
        meas[l] = yrandom[l] + x1[nxx+theta*(l+1)]
        meas[l+L] = yrandom[l+L] + x2[nxx+theta*(l+1)]
        del l
    
    return meas

"""
"""

def measure(r,x,rng,Nx,L,nxx,nyy):
    "function maps the observations into observation space"
    yrandom = np.zeros(int(Nx/10))
    meas = np.zeros(int(Nx/10))
    for i in range(int(Nx/10)):
        meas[i] = yrandom[i] + x[(i+1)*10]
        del i

    return meas

"""

def perturbed_ob(y_t,rng,nens,L,r):
    "generate the perturbed observations"
    y_pert = np.zeros((L,nens))
    for e in range(nens):
        y_pert[:,e] = r*rng.normal(0,1,L)+y_t
        del e

    return y_pert
