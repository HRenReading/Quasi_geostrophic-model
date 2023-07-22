# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 20:03:25 2022

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
from observation import *
from parameter import *


def readStream(filenumber):
    "function to read the unformatted fortran files"
    #read file
    ft = FortranFile(filenumber,'r')
    #put data into a vector
    x = ft.read_reals()
    return x


def real_vec(t_start,t_end,nxx,nyy,nl):
    "read all files of the true state of the system"
    xr = np.zeros((nxx*nyy*nl,t_end-t_start+1))
    for i in range(t_start,t_end+1):
        if i <10: 
            xr[:,i-t_start] = readStream("stream999.0000"+str(i))
        if 10 <= i < 100:
            xr[:,i-t_start] = readStream("stream999.000"+str(i))
        if i >= 100:
            xr[:,i-t_start] = readStream("stream999.00"+str(i))
        del i
    return xr


def read_ensemble(nens,t_start,t_end,nxx,nyy,nl):
    "function read all the ensemble data"
    xb = np.empty((nxx*nyy*nl,nens,t_end-t_start+1)); xb.fill(np.NaN)
    for i in range(1,nens+1):
        for j in range(t_start,t_end+1):
            if i<10:   
                if j<10:   
                    xb[:,i-1,j-t_start] = readStream("stream00"+str(i)+".0000"+str(j))
                if 10<=j<100:
                    xb[:,i-1,j-t_start] = readStream("stream00"+str(i)+".000"+str(j))
                if 100<=j:
                    xb[:,i-1,j-t_start] = readStream("stream00"+str(i)+".00"+str(j))
            if 10<= i<100:
                if j<10:   
                    xb[:,i-1,j-t_start] = readStream("stream0"+str(i)+".0000"+str(j))
                if 10<=j<100:
                    xb[:,i-1,j-t_start] = readStream("stream0"+str(i)+".000"+str(j))
                if 100<=j:
                    xb[:,i-1,j-t_start] = readStream("stream0"+str(i)+".00"+str(j))
            if i>=100:
                if j<10:   
                    xb[:,i-1,j-t_start] = readStream("stream"+str(i)+".0000"+str(j))
                if 10<=j<100:
                    xb[:,i-1,j-t_start] = readStream("stream"+str(i)+".000"+str(j))
                if 100<=j:
                    xb[:,i-1,j-t_start] = readStream("stream"+str(i)+".00"+str(j))
            del j
        del i
    return xb


def vec2field(x,nxx,nyy):
    "function that transfers the vector to field"
    s1 = np.reshape(x[0:nxx*nyy],(nyy,nxx))
    s2 = np.reshape(x[nxx*nyy:2*nxx*nyy],(nyy,nxx))
    s1 = s1*f0/g
    
    return s1,s2 


def real_field(xr,time_step,nxx,nyy,g,f0):
    "generate the field of the real stream function"
    s1_r = np.empty((time_step+1,nyy,nxx)); s1_r.fill(np.NaN)
    s2_r = np.empty((time_step+1,nyy,nxx)); s2_r.fill(np.NaN)
    for i in range(0,time_step+1):
        s1_r[i,:,:],s2_r[i,:,:] = vec2field(xr[:,i],nxx,nyy)
        del i
    return s1_r,s2_r


def en_vec2field(x_en,t,nens,nxx,nyy,g,f0):
    "function transfer the ensemble vector to field"
    s1_en = np.empty((nens,t+1,nyy,nxx)); s1_en.fill(np.NaN)
    s2_en = np.empty((nens,t+1,nyy,nxx)); s2_en.fill(np.NaN)
    for i in range(nens):
        for j in range(t+1):
            s1_en[i,j,:,:],s2_en[i,j,:,:] = vec2field(x_en[:,i,j],nxx,nyy)
            del j
        del i
    return s1_en,s2_en
"""
s1_r,s2_r = real_field(x_r,t,nxx,nyy,g,f0)
s1_en,s2_en = en_vec2field(x_en,t,nens,nxx,nyy,g,f0)     
#s1_a,s2_a = en_vec2field(x_an,t,nens,nxx,nyy,g,f0)  
s1b_mean,s2b_mean = mean_field(s1_en,s2_en,nyy,nxx,nens,t)
#s1a_mean,s2a_mean = mean_field(s1_a,s2_a,nyy,nxx,nens,t)

plt.figure()
plt.subplot(2,1,1)
plt.contourf(s1_r[0,:,:],20,cmap='bwr')
plt.colorbar()
plt.title("Sea-surface height")
plt.show()
plt.subplot(2,1,2)
plt.contourf(s2_r[0,:,:],20,cmap='bwr')
plt.colorbar()
plt.show()


plt.figure()
plt.subplot(3,2,1)
plt.contourf(s1_r[-1,:,:],20,cmap='bwr')
plt.colorbar()
plt.title("Upper layer truth")
plt.show()
plt.subplot(3,2,3)
plt.contourf(s1b_mean[-1,:,:],20,cmap='bwr')
plt.colorbar()
plt.title("Upper layer prior")
plt.show()
plt.subplot(3,2,5)
plt.contourf(s1a_mean[-1,:,:],20,cmap='bwr')
plt.colorbar()
plt.title("Upper layer posterior")
plt.show()
plt.subplot(3,2,2)
plt.contourf(s2_r[-1,:,:],20,cmap='bwr')
plt.colorbar()
plt.title("Lower layer truth")
plt.show()
plt.subplot(3,2,4)
plt.contourf(s2b_mean[-1,:,:],20,cmap='bwr')
plt.colorbar()
plt.title("Lower layer prior")
plt.show()
plt.subplot(3,2,6)
plt.contourf(s2a_mean[-1,:,:],20,cmap='bwr')
plt.colorbar()
plt.title("Lower layer posterior")
plt.show()
"""
def upDimension(x_en,x_r,time,Nx,nens,w_guess):
    'prepare the data for EnKS'
    xb_new = np.zeros((Nx*time+1,nens))
    xr_new = np.zeros(Nx*time)
    for i in range(time):
        xb_new[i*Nx:(i+1)*Nx,:] = x_en[:,:,i]
        xr_new[i*Nx:(i+1)*Nx] = x_r[:,i]
        del i
    xb_new[-1,:] = w_guess
    return xb_new,xr_new


def downDimension(xb,t,Nx,nens):
    'prepare the data for EnKS'
    xb_new = np.zeros((Nx,nens,t+1))
    for i in range(t+1):
        xb_new[:,:,i] = xb[i*Nx:(i+1)*Nx,:]
        del i
    phi_a = xb[-1,:]
    return xb_new,phi_a


def mean_field(s1,s2,nyy,nxx,nens,t):
    "compute the temporal average"
    s1_mean = np.mean(s1,axis=0)
    s2_mean = np.mean(s2,axis=0)
    return s1_mean,s2_mean

def save_data(xr,xb,xa,number_window,nens):
    writeFile('real'+str(number_window),xr)
    for i in range(nens):
        writeFile('prior'+str(number_window)+'.'+str(i),xb[:,i])
        writeFile('post'+str(number_window)+'.'+str(i),xa[:,i])
        del i

