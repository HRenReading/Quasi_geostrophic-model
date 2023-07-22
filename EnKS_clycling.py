# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 13:06:26 2022

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
from scipy.linalg import eigh
from scipy.linalg import circulant
from scipy.stats import gaussian_kde
import math
from math import log
from readData import *
from parameter import *
from observation import *
from RMSE_spread import *
from writeToFiles import *
from IEnKS import *


plt.rc('font', size=18)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=18)    # legend fontsize
plt.rc('figure', titlesize=18)  # fontsize of the figure title

def transfer(phi,nens):
    "convert phi to omega"
    w = np.zeros(nens)
    for e in range(nens):
        w[e] = -1/np.log(phi[e])
        del e
        
    return w

def EnKS(x_r,x_en,t,nxx,nyy,nl,nens,r,rng,L,Nx,phi_g):
    "IEnKS without iteration (EnKS)"
    #obtain the obs
    y_t = measure(r,x_r[:,-1],rng,Nx,L,nxx,nyy)
    #perturbed obs
    y_pert = perturbed_ob(y_t,rng,nens,2*L,r)
    #read data and augmentation
    xb,xr = upDimension(x_en,x_r,t+1,Nx,nens,phi_g) 
    #da process
    xa = IETKS(Nx*(t+1)+1,2*L,nens,xb,y_pert,r,Nx,nxx,nyy,rng,L,x_en[:,:,-1])
    #seperate data and parameters
    x_an,phi_a = downDimension(xa,t,Nx,nens) 
    return x_an,phi_a

def EnKS_cycling(x_r,x_en,t,nxx,nyy,nl,nens,r,rng,L,Nx,phi_g):
    y_t = measure(r,x_r[:,-1],rng,Nx,L,nxx,nyy)
    y_pert = perturbed_ob(y_t,rng,nens,2*L,r)
    xb,xr = upDimension(x_en,x_r,t,Nx,nens,phi_g) 
    xa = IETKS(Nx*t+1,2*L,nens,xb,y_pert,r,Nx,nxx,nyy,rng,L,x_en[:,:,-1])
    x_an,phi_a = downDimension(xa,t-1,Nx,nens)
    return x_an,phi_a

def plot_real(xr,t,nxx,nyy,g,f0,nwin_step):
    "plot the real field at the end of every iteration"
    s1,s2 = real_field(xr,t,nxx,nyy,g,f0)
    time = [0,int(t/2),t]
    plt.figure()
    for i in range(3):
        plt.subplot(3,1,i+1)
        plt.contourf(s1[time[i],:,:],20,cmap='bwr')
        plt.colorbar()
        plt.title('Upper layer truth at t = '+str(time[i]+nwin_step))
        plt.show()
        del i
    plt.figure()
    for i in range(3):
        plt.subplot(3,1,i+1)
        plt.contourf(s2[time[i],:,:],20,cmap='bwr')
        plt.colorbar()
        plt.title('Lower layer truth at t = '+str(time[i]+nwin_step))
        plt.show()
        del i
        
def phiTow(phi,nens,n_windows):
    "convert phi to omega"
    w = np.zeros((nens,n_windows))
    for i in range(n_windows):
        for e in range(nens):
            w[e,i] = -1/np.log(phi[e,i])
            del e
        del i
    return w


def data_plot(h):
    #generate the normal/Guassian plot
    h.sort()
    hmean = np.mean(h)
    hstd = np.std(h)
    pdf = stats.norm.pdf(h, hmean, hstd)
    return pdf

        
def plot_pdf(w_guess,phi_a,nens):
    w_a = transfer(phi_a,nens)
    plt.figure()
    plt.xlabel('$\omega$')
    plt.plot(w_guess,data_plot(w_guess),color = 'b',label = '$\omega^b$')
    plt.plot(w_a,data_plot(w_a),color = 'r',label = '$\omega^a$')
    plt.axvline(1.0,color = 'black',label = '$\omega_r$')
    plt.xlabel('$\omega$')
    plt.ylabel('Probability density')
    plt.legend(loc = 'upper right')
    plt.show()
    
def select(phi,w,nens):
    #select useful parameter
    for i in range(nens):
        if phi[i] < 0:
            phi[i] = phi[i-1]
        else:
            phi[i] = phi[i]
        if w[i] >= 5:
            w[i] = w[i-1]
        else:
            w[i] = w[i]
        del i
    return phi,w


#matrix to store parameter
w_all = np.zeros((nens,n_windows))
#output parameter
writeFile("phi_guess",phi_g)
#call Fortran code to run the QG model in Fortarn (generate data)
os.system("gfortran -o QG -fdefault-real-8  -fdefault-double-8 -fallow-invalid-boz HFFT.f nagsources.f parameter.f90  random.f90 helm.f90 PV.f90 QGModel.f90 -lBLAS -lLAPACK")
os.system("QG")
#read data
x_r = real_vec(0,t,nxx,nyy,nl)     #real state of the system in vector
x_en = read_ensemble(nens,0,t,nxx,nyy,nl) #ensemble members
#perform DA (EnKS)
x_an,phi_a = EnKS(x_r,x_en,t,nxx,nyy,nl,nens,r,rng,L,Nx,phi_g)
ini_r,ini = initial(x_r[:,-1],x_an[:,:,-1],nxx,nyy,nl,nens) #obtain initial condition for the next cycle
#write initials to unformatted file
writeini(ini_r,ini,nens)
#convert the parameter
w_a = transfer(phi_a,nens)
phi_a,w_a = select(phi_a,w_a,nens)
w_all[:,0] = w_a
#write to file
writeFile("phi_guess",phi_a)
#cycling for multiple windows
for i in range(1,n_windows):
    #write the analyzed data at the end of the simulation as the initial
    os.system("gfortran -o QGCycling -fdefault-real-8  -fdefault-double-8 -fallow-invalid-boz HFFT.f nagsources.f parameter.f90  random.f90 helm.f90 PV.f90 QGCycling.f90 -lBLAS -lLAPACK")
    os.system("QGCycling") 
    x_r = real_vec(1,t,nxx,nyy,nl) 
    x_en = read_ensemble(nens,1,t,nxx,nyy,nl)
    x_an,phi_a = EnKS_cycling(x_r,x_en,t,nxx,nyy,nl,nens,r,rng,L,Nx,phi_a)
    ini_r,ini = initial(x_r[:,-1],x_an[:,:,-1],nxx,nyy,nl,nens)
    writeini(ini_r,ini,nens)
    w_a = transfer(phi_a,nens)
    phi_a,w_a = select(phi_a,w_a,nens)
    w_all[:,i] = w_a
    writeFile("phi_guess",phi_a)
    del i 




"""

def readwindows(Nx,n_windows,nens):
    "read saved data"
    xr = np.zeros((Nx,n_windows))
    xb = np.zeros((Nx,nens,n_windows))
    xa = np.zeros((Nx,nens,n_windows))
    for i in range(n_windows):
        xr[:,i] = readStream("real"+str(i))
        for e in range(nens):
            xb[:,e,i] = readStream("prior"+str(i)+'.'+str(e))
            xa[:,e,i] = readStream("post"+str(i)+'.'+str(e))
            del e
        del i
    return xr,xb,xa


"""

colors = plt.cm.Reds(np.linspace(0,1,5))

plt.figure()
plt.xlabel('$\omega$')
plt.plot(w_guess,data_plot(w_guess),color = 'b',label = '$\omega^b$')
for i in range(n_windows):
    plt.plot(w_all[:,i],data_plot(w_all[:,i]),color = colors[i-5],label = '$\omega^a$ at '+str(i+1)+' windows')
    del i
plt.axvline(1.0,color = 'black')
plt.xlabel('$\omega$')
plt.ylabel('Probability density')
plt.legend(loc = 'upper right')
plt.show()

"""

RMSE_1b_all = np.zeros(n_windows*t+1)
spread_1b_all = np.zeros(n_windows*t+1)
RMSE_1a_all = np.zeros(n_windows*t+1)
spread_1a_all = np.zeros(n_windows*t+1)
RMSE_2b_all = np.zeros(n_windows*t+1)
spread_2b_all = np.zeros(n_windows*t+1)
RMSE_2a_all = np.zeros(n_windows*t+1)
spread_2a_all = np.zeros(n_windows*t+1)

rmse1_b,rmse2_b = RMSE(x_en,x_r,t+1,nyy,nxx,nens,f0,g)
rmse1_a,rmse2_a = RMSE(x_an,x_r,t+1,nyy,nxx,nens,f0,g)
spread1_a,spread2_a = spread(x_an,t+1,nyy,nxx,nens,f0,g)
spread1_b,spread2_b = spread(x_en,t+1,nyy,nxx,nens,f0,g)

RMSE_1b_all[0:t+1] = rmse1_b
spread_1b_all[0:t+1] = spread1_b
RMSE_1a_all[0:t+1] = rmse1_a
spread_1a_all[0:t+1] = spread1_a
RMSE_2b_all[0:t+1] = rmse2_b
spread_2b_all[0:t+1] = spread2_b
RMSE_2a_all[0:t+1] = rmse2_a
spread_2a_all[0:t+1] = spread2_a

n_windows=1
plt.figure()
plt.subplot(2,2,1)
plt.plot(np.arange(0,t*n_windows+1),rmse1_a,color = 'r',label = 'Posterior RMSE')
plt.plot(np.arange(0,t*n_windows+1),rmse1_b,color = 'b',label = 'Prior RMSE')
plt.title('Upper layer RMSE')
plt.show()
plt.legend(loc='upper left')
plt.subplot(2,2,2)
plt.plot(np.arange(0,t*n_windows+1),rmse2_a,color = 'r',label = 'Posterior RMSE')
plt.plot(np.arange(0,t*n_windows+1),rmse2_b,color = 'b',label = 'Prior RMSE')
plt.title('Lower layer RMSE')
plt.show()
plt.legend(loc='upper left')
plt.subplot(2,2,3)
plt.plot(np.arange(0,t*n_windows+1),spread1_a,color = 'r',label = 'Posterior spread')
plt.plot(np.arange(0,t*n_windows+1),spread1_b,color = 'b',label = 'Prior spread')
plt.title('Upper layer spread')
plt.show()
plt.legend(loc='upper left')
plt.subplot(2,2,4)
plt.plot(np.arange(0,t*n_windows+1),spread2_a,color = 'r',label = 'Posterior spread')
plt.plot(np.arange(0,t*n_windows+1),spread2_b,color = 'b',label = 'Prior spread')
plt.title('Lower layer spread')
plt.show()
plt.legend(loc='upper left')
"""

