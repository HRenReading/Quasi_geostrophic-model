# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 17:55:27 2022

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
import scipy.stats as stats
from scipy.linalg import circulant
from scipy.stats import gaussian_kde
import math
from readData import *
from parameter import *
from observation import *
from RMSE_spread import *


def writeFile(filename,variable):
    """
    f = FortranFile(filename, mode='w')
    f.write_record(variable)
    f.close()
    """
    f = FortranFile(filename, mode='w')
    data = np.transpose(variable)
    f.write_record(data)
    f.close()


def writeini(ini_r,ini_b,nens):
    'write data to Fortran unformatted'
    writeFile('ini_r',ini_r)
    for i in range(1,nens+1):
        if i < 10:
            writeFile('ini00'+str(i),ini_b[:,:,:,i-1])
        if 10 <= i < 100:
            writeFile('ini0'+str(i),ini_b[:,:,:,i-1])
        if i >= 100:
            writeFile('ini'+str(i),ini_b[:,:,:,i-1])

    
def initial(xr,xa,nxx,nyy,nl,nens):
    ini = np.zeros((nxx,nyy,nl,nens))
    ini_r = np.zeros((nxx,nyy,nl))
    s1_r,s2_r = vec2field(xr,nxx,nyy) 
    s1_a = np.zeros((nens,nyy,nxx))
    s2_a = np.zeros((nens,nyy,nxx))
    for e in range(nens):
        s1_a[e,:,:],s2_a[e,:,:] = vec2field(xa[:,e],nxx,nyy)
        del e
    ini_r[:,:,0] = s1_r[:,:].T
    ini_r[:,:,1] = s2_r[:,:].T
    for i in range(nens):
        ini[:,:,0,i] = s1_a[i,:,:].T
        ini[:,:,1,i] = s2_a[i,:,:].T
        del i
    return ini_r,ini
