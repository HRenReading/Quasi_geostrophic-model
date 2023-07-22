# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 13:09:02 2022

@author: 44754
"""

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


def RMSE(xb,xr,t,nyy,nxx,nens,f0,g):
    "compute the ensemble RMSE"
    l1 = xb[0:nxx*nyy,:,:]*f0/g
    l2 = xb[nxx*nyy:,:,:]
    xr1 = xr[0:nxx*nyy,:]*f0/g
    xr2 = xr[nxx*nyy:,:]
    l1_mean = np.zeros((nxx*nyy,t))
    l2_mean = np.zeros((nxx*nyy,t))
    for i in range(nens):
        l1_mean += l1[:,i,:]
        l2_mean += l2[:,i,:]
        del i
    l1_mean = l1_mean/nens
    l2_mean = l2_mean/nens
    rmse1 = np.zeros(t)
    rmse2 = np.zeros(t)
    for i in range(t):
        rmse1[i] = np.sqrt(np.dot((xr1[:,i]-l1_mean[:,i]).T,(xr1[:,i]-l1_mean[:,i]))/(nxx*nyy))
        rmse2[i] = np.sqrt(np.dot((xr2[:,i]-l2_mean[:,i]).T,(xr2[:,i]-l2_mean[:,i]))/(nxx*nyy))
        del i
       
    return rmse1,rmse2


def spread(xb,t,nyy,nxx,nens,f0,g):
    "compute the ensemble RMSE"
    l1 = xb[0:nxx*nyy,:,:]*f0/g
    l2 = xb[nxx*nyy:,:,:]
    l1_mean = np.zeros((nxx*nyy,t))
    l2_mean = np.zeros((nxx*nyy,t))
    for i in range(nens):
        l1_mean += l1[:,i,:]
        l2_mean += l2[:,i,:]
        del i
    l1_mean = l1_mean/nens
    l2_mean = l2_mean/nens
    spread1 = np.zeros(t)
    spread2 = np.zeros(t)
    sum1 = np.zeros(t)
    sum2 = np.zeros(t)
    for i in range(t):
        for j in range(nens):
            sum1[i] += (np.dot((l1[:,j,i]-l1_mean[:,i]).T,(l1[:,j,i]-l1_mean[:,i]))/(nens-1))/(nyy*nxx)
            sum2[i] += (np.dot((l2[:,j,i]-l2_mean[:,i]).T,(l2[:,j,i]-l2_mean[:,i]))/(nens-1))/(nyy*nxx)
            del j
        del i
   
    for i in range(t):
        spread1[i] = np.sqrt(sum1[i])
        spread2[i] = np.sqrt(sum2[i])
        del i
     
    return spread1,spread2