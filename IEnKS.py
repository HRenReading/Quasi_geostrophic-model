# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 17:56:57 2022

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
from scipy.linalg import eigh
import math
from observation import *


def IETKS(nx,ny,nens,xp,y_pert,sdobs,Nx,nxx,nyy,rng,L,x_en): 
    "EnKS in subspace with iteration as a option"
    gamma = 1.0                     # iteration step size
    ymodel = np.zeros((ny,nens))    # Observation forecast ensemble perturbations
    w = np.zeros((nens,nens))
    scale = 1./np.sqrt(nens-1.)     # scale factor for ensemble matrices  
    Imean = (np.identity(nens) - np.ones((nens,nens))/nens)*scale
    for e in range(nens):
        ymodel[:,e] = measure(0,x_en[:,e],rng,Nx,L,nxx,nyy)
        del e
    ymodel_pert = np.dot(ymodel,Imean)
    #ymodel_pert = np.dot(ymodel_pert,np.dot(np.linalg.pinv(AA),AA))
    Momega = np.identity(nens) + np.dot(w,Imean)
    s = np.dot(ymodel_pert,np.linalg.pinv(Momega))
    D_tilt = np.dot(s,w) + y_pert - ymodel
    u, s, vh = np.linalg.svd(s, full_matrices=False)
    s2 = np.zeros((nens,nens))
    s2[0,0] = 1./s[0]
    tol = 1.e-3
    for i in range(1,nens):
        if s[i] > tol * s[1]:
          s2[i,i] = 1./s[i]
        else:
          s2[i,i] = 0

    a = np.dot(s2, np.dot(u.T,ymodel_pert))
    a =np.dot(a,a.T)
    l,Q = eigh(a) 
    l = np.diag(l)
    p1 = np.dot(vh.T,Q)
    p2 = np.linalg.pinv(np.identity(nens)+l)
    p3 = np.dot(Q.T,np.dot(s2,u.T))
    P1 = np.dot(p1,p2)
    P2 = np.dot(P1,p3)
    P3 = np.dot(P2,D_tilt)
    w = w - gamma*(w - P3)
    T = np.identity(nens) + w*scale
    xp = np.dot(xp,T)


    return xp

"""

def IETKS(nx,ny,nens,xp,y_pert,sdobs,Nx,nxx,nyy,rng,L,x_en):
    "ensemble subspace RML"
    gamma = 1.0
    W = np.zeros((nens,nens))
    scale = 1./np.sqrt(nens-1.)
    Imean = (np.identity(nens) - np.ones((nens,nens))/nens)*scale
    E = np.dot(y_pert,Imean)
    ymodel = np.zeros((ny,nens))
    for e in range(nens):
        ymodel[:,e] = measure(0,x_en[:,e],rng,Nx,L,nxx,nyy)
        del e
    ymodel_hat = np.dot(ymodel,Imean)
    Omega = np.identity(nens) + np.dot(W,Imean)
    S = np.dot(ymodel_hat,np.linalg.pinv(Omega))
    D_hat = np.dot(S,W)+y_pert-ymodel
    inva = np.dot(S,S.T)+np.dot(E,E.T)
    w = W - gamma*(W - np.dot(np.dot(S.T,np.linalg.pinv(inva)),D_hat))
    T = np.identity(nens)+w*scale
    xa = np.dot(xp,T)
    
    return xa


def IETKS(nx,ny,nens,xp,y_pert,sdobs,Nx,nxx,nyy,rng,L,x_en):    
    ymodel = np.zeros([ny,nens])    # Observation forecast ensemble perturbations
    H      = np.zeros([ny,nens])    # Intermediate matrix
    s      = np.zeros([ny,nens])    # Intermediate matrix

    IIy    = np.zeros([ny,ny])     # identity matrix in observation space
    a      = np.zeros([ny,ny])     # Intermediate matrix
   
    w     = np.zeros([nens,nens])    # ensemble weight matrix
    Imean = np.zeros([nens,nens])    # multiply to obtain ensemble mean
    IIw   = np.zeros([nens,nens])    # identity matrix in ensemble space
    Momega = np.zeros([nens,nens])    # Intermediate matrix
    T     = np.zeros([nens,nens])    # Intermediate matrix
    AA    = np.copy(xp)
   
    scale = 1./np.sqrt(nens-1.)      # scale factor for ensemble matrices  
    gamma = 1.0                   # iteration step size
    IIy = np.identity(ny)
    IIw = np.identity(nens)
    Imean = (IIw - np.ones((nens,nens))/nens)*scale
    for e in range(nens):
        ymodel[:,e] = measure(0,x_en[:,e],rng,Nx,L,nxx,nyy)
        del e
    ymodel_pert = np.dot(ymodel,Imean)
    AA = np.dot(xp,Imean)
    #ymodel_pert = np.dot(ymodel_pert,np.dot(np.linalg.pinv(AA),AA))
    Momega = IIw + np.dot(w,Imean)
    s = np.dot(ymodel_pert,np.linalg.inv(Momega))
    H = np.dot(s,w)+ y_pert - ymodel
    a = np.dot(s,s.T)+sdobs*sdobs*IIy
    w = w - gamma*(w - np.dot(np.transpose(s),np.dot(np.linalg.pinv(a),H)))
    T = (IIw + w*scale)
    xp = np.dot(xp,T)

    return xp
"""
