# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 20:05:09 2022

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

random.seed(1)
np.random.seed(1) 
rng = rnd.default_rng(seed=100)


#general parameters
nxx = 257              #zonal (x) dimension of the system
nyy = 129              #medridional (y) dimension of the system
nl = 2                 #number of vertical layers 
nens = 500             #number of ensemble members
Nx = nxx*nyy*nl        #number of points of the model
n_it = 10              #number of iterations
t = 10           #number of time steps in one simulation window
g = 9.81       #gravity acceleration
f0 = -1.26e-4  #Coriolis coefficient
w_r = 1.0
w_g = 3.0
w_guess = abs(np.random.normal(w_g,1.,nens))
phi_g = np.exp(-1/w_guess)
n_windows = 5
theta = 5
L = int(nxx-4)*int(nyy)
r = 5.e3
