#!/usr/bin/env python
# -*- coding: utf-8 -*-

# $Id: sw_1_4_assim.py 520 2015-11-23 09:18:24Z michael.wuersch $

# based upon script from Yvonne and alterated by Kirsten

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pdb
from matplotlib import pyplot as plt
import math
from constants import *
from Plotting import rmse
from DA_2019 import *
from scipy import stats
from datetime import datetime
from numpy import ndarray
from netCDF4 import Dataset
import numpy.ma as ma
import pickle
from datetime import datetime
from scipy.stats import norm

# note those functions appended with a '__' use every timestep. 
# msw_model_expanded_3 includes domain dimension of phi_c. 
# msw_model_expanded_4 includes new initialisation seet up. DA is 1 24hour diurnal cycle. The 40 members are all the same at beginning of DA. 

def sweq__(x,phi_data,phi_c_data,nsteps,n,n_ens,mindex,nindex,time,batch_number,members_per_batch):
    """
    x is the concatenated state of size 3*n. This is the relevant model input and ouput.
    nsteps is the number of timesteps
    n is the total grid size
    n_ens is the number of ensemble members
    mindex = mountain index - 1: bell shaped mountain
    nindex = noise index - 1: stochastic Gaussian noise throughout grid domain
    time indicates whether on first time step. Used to install orography and keep count for random seed
    """  
    phi = np.zeros((n+2,n_ens))                                 # geopotential
    phi_c = np.zeros((n,3,n_ens))                                 # constant geopotential when h>h_c
    beta_new = np.zeros((n+2,n_ens))                            # returns matrix of n+2 columns. Shape (n+2)
    harray = np.zeros((n+2,n_ens))                              # array which will be able to alter the depth of fluid
    y = np.zeros((3*n,n_ens,nsteps))                            # final outputting array
    y_phi_c = np.zeros((n,2,n_ens,nsteps))                        # final outputting array
    y_phi = np.zeros((n,n_ens,nsteps))                          # final outputting array
    array = np.zeros((n,n_ens))                                 # used in calculations

    ku = 2000                                                   # diffusion constant for u usually 2000
    kh = 6000                                                   # diffusion coefficient for depth h usually 6000
    
    if mindex == 1:                                             # bell shaped mountain ridge
      ku = 3000                                                 # diffusion constant for u
      kh = ku                                                   # diffusion coefficient for depth h
      for i in range(n):
        harray[i+1,:] = amp*hw*hw/((float(i)-mu)*(float(i)-mu)+hw*hw)   # every ensemble member has same mountain

    if mindex == 2:                                             # k**(-1) power spectrum
      x_sin = np.arange(0,n)                                    # mindex #1 and #2 need work and checking
      s = 0.2+np.sin((x_sin-100)/500)
      np.random.seed(2)
      data = np.random.rand(n)-0.5                              # data centred around zero
      ps = np.abs(np.fft.fft(data,norm='ortho'))**2             # calculating the power spectrum
      for i in range(n_ens):
        harray[1:n+1,i]=ps*s

    harray[0,:] = harray[n,:]
    harray[n+1,:] = harray[1,:]

    u = np.zeros((3,n+2,n_ens))                                             # horizontal velocity returns matrix of 3 row and n+2 columns 
    h = np.zeros((3,n+2,n_ens))                                             # depth, returns matrix of 3 row and n+2 columns 
    r = np.zeros((3,n+2,n_ens))                                             # rain, returns matrix of 3 row and n+2 columns 
    u[0,1:n+1,:],u[1,1:n+1,:] = x[0:n,:],x[0:n,:]                           # filling in the middle of the u,h and r matrixes with input 
    h[0,1:n+1,:],h[1,1:n+1,:] = x[n:2*n,:],x[n:2*n,:]
    r[0,1:n+1,:],r[1,1:n+1,:] = x[2*n:3*n,:],x[2*n:3*n,:]
    
    u[1,0,:],u[0,0,:] = u[0,n,:],u[0,n,:]                                   # filling in the 2 outer columns of the u, h and r matrixes
    u[1,n+1,:],u[0,n+1,:] = u[0,1,:],u[0,1,:]                               # for boundary conditions from which the model effectively works. 
    h[1,0,:],h[0,0,:] = h[0,n,:],h[0,n,:]
    h[1,n+1,:],h[0,n+1,:] =  h[0,1,:],h[0,1,:]
    r[1,0,:],r[0,0,:] = r[0,n,:],r[0,n,:]
    r[1,n+1,:],r[0,n+1,:] = r[0,1,:],r[0,1,:]

    phi[1:n+1,:] = phi_data   # now fill the phi arrays with that data
    phi_c[:,:2,:] = phi_c_data   # now fill the phi arrays with that data
    
    phi[0,:] = phi[n,:] # do periodicness as for u,h and r above
    phi[n+1,:] = phi[1,:]

    if time == 0:                                                           # intall ridge - only happens at first time step
      h[:,:,:] = h[:,:,:] - harray[:,:]     

    for it in range(nsteps):                                                # loop over the given model timestep. Saves output after nsteps
      if nindex == 1:
        noise_seed = time+it
        u = noise(u,n,n_ens,n_array,noise_seed,batch_number,members_per_batch)                             # returns the u matrix with the 2nd row with stochastic noise across the ensemble
        u[1,0,:] = u[1,n,:]
        u[1,n+1,:] = u[1,1,:]
      
      array = np.where(h[1,1:n+1,:]<90.4,0,(u[1,2:n+2,:] - u[1,1:n+1,:])/dx)
      array = np.where(array[:,:]>0,0,array[:,:])
      con_av = np.mean(array[:,:],axis=0)
      con_av = np.nan_to_num(con_av)                   

      S_conv = gamma_2 * beta * con_av[:]
     
      phi_c[:,2,:] = phi_c[:,0,:] + 2*dts*(S_rad + S_for + S_conv)

      if mindex == 1 or mindex==2:
        phi[1:n+1,:] = np.where( h[1,1:n+1,:]+harray[1:n+1,:] > h_cloud , phi_c[:,2,:]+g*harray[1:n+1,:], g*h[1,1:n+1,:] ) # if condition met, return phic in h matrix. If not, return g*h thing which would be the normal geopotential below the first threshold
      else:
        phi[1:n+1,:] = np.where( h[1,1:n+1,:] > h_cloud , phi_c[:,2,:], g*h[1,1:n+1,:] )

      phi[0,:]   = phi[n,:]                                             # boundary conditions for geopotential
      phi[n+1,:] = phi[1,:]
      phi[:,:] = phi[:,:] + gamma * r[1,:,:] 

      # shallow water equations =D
      u[2,1:n+1,:] = u[0,1:n+1,:] - (dts/(2*dx))*(u[1,2:n+2,:]**2 - u[1,0:n,:]**2) - (2*dts/dx)*(phi[1:n+1,:]-phi[0:n,:])  + (ku/(dx*dx))*(u[0,2:n+2,:] - 2*u[0,1:n+1,:] + u[0,0:n,:])*dts*2     # momentum equation  # fixed diffusion term
      h[2,1:n+1,:] = h[0,1:n+1,:] - (dts/dx)*(u[1,2:n+2,:]*(h[1,1:n+1,:]+h[1,2:n+2,:]) - u[1,1:n+1,:]*(h[1,0:n,:]+h[1,1:n+1,:])) + (kh/(dx*dx))*(h[0,2:n+2,:] - 2*h[0,1:n+1,:] + h[0,0:n,:])*dts*2   # continuity equation  # fixed diffusion term

      if mindex == 1 or mindex==2:
        mask = np.logical_and(h[1,1:n+1,:]+harray[1:n+1,:] > h_rain, u[1,2:n+2,:]-u[1,1:n+1,:] < 0)    # conditions for rain
      else:
        mask = np.logical_and(h[1,1:n+1,:] > h_rain, u[1,2:n+2,:]-u[1,1:n+1,:] < 0)   

      beta_new[1:n+1,:] = np.where( mask, beta , 0 )
      r[2,1:n+1,:] = r[0,1:n+1,:] - (dts/(2*dx))*(u[1,2:n+2,:]+u[1,1:n+1,:])*(r[1,2:n+2,:]-r[1,0:n,:]) - alpha*dts*2.0*r[1,1:n+1,:]-2.0*beta_new[1:n+1,:]*(dts/dx)*(u[1,2:n+2,:]-u[1,1:n+1,:]) + (kr/(dx*dx))*(r[0,2:n+2,:] - 2.0*r[0,1:n+1,:] + r[0,0:n,:])*dts*2 # rain equation  # with advection # always peaky
      #r[2,1:n+1,:] = r[0,1:n+1,:] - alpha*dts*2.0*r[1,1:n+1,:]-2.0*beta_new[1:n+1,:]*(dts/dx)*(u[1,2:n+2,:]-u[1,1:n+1,:]) + (kr/(dx*dx))*(r[0,2:n+2,:] - 2.0*r[0,1:n+1,:] + r[0,0:n,:])*dts*2 # rain equation, no advection
      #r[2,1:n+1,:] = r[0,1:n+1,:] - alpha*dts*2.0*r[1,1:n+1,:]-2.0*beta_new[1:n+1,:]*(dts/dx)*(u[1,2:n+2,:]-u[1,1:n+1,:]) + (kr/(dx*dx))*(r[0,2:n+2,:] - 2.0*r[0,1:n+1,:] + r[0,0:n,:])*dts*2  - (dts/(dx))*(u[1,1:n+1,:])*(r[1,2:n+2,:]-r[1,0:n,:])  # other advection term 
     
      u[2,0,:]   = u[2,n,:]                                                  # boundary conditions
      u[2,n+1,:] = u[2,1,:]
      h[2,0,:]   = h[2,n,:]
      h[2,n+1,:] = h[2,1,:]
      r[2,0,:]   = r[2,n,:]
      r[2,n+1,:] = r[2,1,:]
      r[2,:,:] = np.where(r[2,:,:]<r_threshold,0,r[2,:,:])

      d = filter_parameter*.5*(u[2,:,:] - 2.*u[1,:,:] + u[0,:,:])            # RAW filter. Accounts for the growing computational mode.
      u[0,:,:] = u[1,:,:] + alpha_filt*d
      u[1,:,:] = u[2,:,:] - (1-alpha_filt)*d                     
      d = filter_parameter*.5*(h[2,:,:] - 2.*h[1,:,:] + h[0,:,:])
      h[0,:,:] = h[1,:,:] + alpha_filt*d
      h[1,:,:] = h[2,:,:] - (1-alpha_filt)*d
      d = filter_parameter*.5*(r[2,:,:] - 2.*r[1,:,:] + r[0,:,:])
      r[0,:,:] = r[1,:,:] + alpha_filt*d
      r[1,:,:] = r[2,:,:] - (1-alpha_filt)*d

      d = filter_parameter*.5*(phi_c[:,2,:] - 2.*phi_c[:,1,:] + phi_c[:,0,:])      # is this required here for phi_c???
      phi_c[:,0,:] = phi_c[:,1,:] + alpha_filt*d
      phi_c[:,1,:] = phi_c[:,2,:] - (1-alpha_filt)*d

      y[0:n,:,it] = u[2,1:n+1,:]
      y[n:2*n,:,it] = h[2,1:n+1,:]
      y[2*n:3*n,:,it] = r[2,1:n+1,:] 
      y_phi_c[:,:,:,it] =  phi_c[:,[1,2],:]
      y_phi[0:n,:,it] =  phi[1:n+1,:]

    return y, y_phi_c, y_phi

def sweqss__(x,phi_data,phi_c_data,nsteps,n,n_ens,mindex,nindex,time,batch_number,members_per_batch): # ss for selected saving 
    """
    x is the concatenated state of size 3*n. This is the relevant model input and ouput.
    nsteps is the number of timesteps
    n is the total grid size
    n_ens is the number of ensemble members
    mindex = mountain index - 1: bell shaped mountain
    nindex = noise index - 1: stochastic Gaussian noise throughout grid domain
    time indicates whether on first time step. Used to install orography and keep count for random seed
    """  
    phi = np.zeros((n+2,n_ens))                                 # geopotential
    phi_c = np.zeros((n,3,n_ens))                                 # constant geopotential to use when h>h_c
    beta_new = np.zeros((n+2,n_ens))                            # returns matrix of n+2 columns. Shape (n+2)
    harray = np.zeros((n+2,n_ens))                              # array which will be able to alter the depth of fluid
    y = np.zeros((3*n,n_ens))                                   # final outputting array
    y_phi_c = np.zeros((n,2,n_ens))                               # final outputting arrays
    y_phi = np.zeros((n,n_ens))
    array = np.zeros((n,n_ens))                                 # to be used in calculation

    ku = 2000                                                   # diffusion constant for u, usually 2000
    kh = 6000                                                   # diffusion coefficient for depth h, usually 6000
    
    if mindex == 1:                                             # bell shaped mountain ridge
      ku = 3000                                                 # diffusion constant for u
      kh = ku                                                   # diffusion coefficient for depth h
      for i in range(n):
        harray[i+1,:] = amp*hw*hw/((float(i)-mu)*(float(i)-mu)+hw*hw)   # every ensemble member has same mountain

    if mindex == 2:                                             # k**(-1) power spectrum
      x_sin = np.arange(0,n)                                    # mindex #1 and #2 need work and checking
      s = 0.2+np.sin((x_sin-100)/500)
      np.random.seed(2)
      data = np.random.rand(n)-0.5                              # data centred around zero
      ps = np.abs(np.fft.fft(data,norm='ortho'))**2             # calculating the power spectrum
      for i in range(n_ens):
        harray[1:n+1,i]=ps*s

    harray[0,:] = harray[n,:]
    harray[n+1,:] = harray[1,:]

    u = np.zeros((3,n+2,n_ens))                                             # horizontal velocity returns matrix of 3 row and n+2 columns 
    h = np.zeros((3,n+2,n_ens))                                             # depth, returns matrix of 3 row and n+2 columns 
    r = np.zeros((3,n+2,n_ens))                                             # rain, returns matrix of 3 row and n+2 columns 
    u[0,1:n+1,:],u[1,1:n+1,:] = x[0:n,:],x[0:n,:]                           # filling in the middle of the u,h and r matrixes with input 
    h[0,1:n+1,:],h[1,1:n+1,:] = x[n:2*n,:],x[n:2*n,:]
    r[0,1:n+1,:],r[1,1:n+1,:] = x[2*n:3*n,:],x[2*n:3*n,:]
    
    u[1,0,:],u[0,0,:] = u[0,n,:],u[0,n,:]                                   # filling in the 2 outer columns of the u, h and r matrixes
    u[1,n+1,:],u[0,n+1,:] = u[0,1,:],u[0,1,:]                               # for boundary conditions from which the model effectively works. 
    h[1,0,:],h[0,0,:] = h[0,n,:],h[0,n,:]
    h[1,n+1,:],h[0,n+1,:] =  h[0,1,:],h[0,1,:]
    r[1,0,:],r[0,0,:] = r[0,n,:],r[0,n,:]
    r[1,n+1,:],r[0,n+1,:] = r[0,1,:],r[0,1,:]

    phi[1:n+1,:] = phi_data
    phi_c[n,:2,:] = phi_c_data
    
    phi[0,:] = phi[n,:] # do periodicness as with u,h and r above
    phi[n+1,:] = phi[1,:]

    if time == 0:                                                           # intall ridge - only happens at first time step
      h[:,:,:] = h[:,:,:] - harray[:,:]     

    for it in range(nsteps):                                                # loop over the given model timestep. Saves output after nsteps
      if nindex == 1:
        noise_seed = time+it
        u = noise(u,n,n_ens,n_array,noise_seed,batch_number,members_per_batch)                             # returns the u matrix with the 2nd row with stochastic noise across the ensemble
        u[1,0,:] = u[1,n,:]
        u[1,n+1,:] = u[1,1,:]

      array = np.where(h[1,1:n+1,:]<90.4,0,(u[1,2:n+2,:] - u[1,1:n+1,:])/dx)
      array = np.where(array[:,:]>0,0,array[:,:])
      con_av = np.mean(array[:,:],axis=0)
      con_av = np.nan_to_num(con_av)

      S_conv = gamma_2 * beta * con_av[:]

      phi_c[:,2,:] = phi_c[:,0,:] + 2*dts*(S_rad + S_for + S_conv)

      if mindex == 1 or mindex==2:
        phi[1:n+1,:] = np.where( h[1,1:n+1,:]+harray[1:n+1,:] > h_cloud , phi_c[:,2,:]+g*harray[1:n+1,:], g*h[1,1:n+1,:] ) # if condition met, return phic in h matrix. If not, return g*h thing which would be the normal geopotential below the first threshold
      else:
        phi[1:n+1,:] = np.where( h[1,1:n+1,:] > h_cloud , phi_c[:,2,:], g*h[1,1:n+1,:] )

      phi[0,:]   = phi[n,:]                                             # boundary conditions for geopotential
      phi[n+1,:] = phi[1,:]
      phi[:,:] = phi[:,:] + gamma * r[1,:,:] 

      # shallow water equations =D
      u[2,1:n+1,:] = u[0,1:n+1,:] - (dts/(2*dx))*(u[1,2:n+2,:]**2 - u[1,0:n,:]**2) - (2*dts/dx)*(phi[1:n+1,:]-phi[0:n,:])  + (ku/(dx*dx))*(u[0,2:n+2,:] - 2*u[0,1:n+1,:] + u[0,0:n,:])*dts*2     # momentum equation  # fixed diffusion term
      h[2,1:n+1,:] = h[0,1:n+1,:] - (dts/dx)*(u[1,2:n+2,:]*(h[1,1:n+1,:]+h[1,2:n+2,:]) - u[1,1:n+1,:]*(h[1,0:n,:]+h[1,1:n+1,:])) + (kh/(dx*dx))*(h[0,2:n+2,:] - 2*h[0,1:n+1,:] + h[0,0:n,:])*dts*2   # continuity equation  # fixed diffusion term

      if mindex == 1 or mindex==2:
        mask = np.logical_and(h[1,1:n+1,:]+harray[1:n+1,:] > h_rain, u[1,2:n+2,:]-u[1,1:n+1,:] < 0)    # conditions for rain
      else:
        mask = np.logical_and(h[1,1:n+1,:] > h_rain, u[1,2:n+2,:]-u[1,1:n+1,:] < 0)   

      beta_new[1:n+1,:] = np.where( mask, beta , 0 )
      r[2,1:n+1,:] = r[0,1:n+1,:] - (dts/(2*dx))*(u[1,2:n+2,:]+u[1,1:n+1,:])*(r[1,2:n+2,:]-r[1,0:n,:]) - alpha*dts*2.0*r[1,1:n+1,:]-2.0*beta_new[1:n+1,:]*(dts/dx)*(u[1,2:n+2,:]-u[1,1:n+1,:]) + (kr/(dx*dx))*(r[0,2:n+2,:] - 2.0*r[0,1:n+1,:] + r[0,0:n,:])*dts*2 # rain equation  # with advection # always peaky
      #r[2,1:n+1,:] = r[0,1:n+1,:] - alpha*dts*2.0*r[1,1:n+1,:]-2.0*beta_new[1:n+1,:]*(dts/dx)*(u[1,2:n+2,:]-u[1,1:n+1,:]) + (kr/(dx*dx))*(r[0,2:n+2,:] - 2.0*r[0,1:n+1,:] + r[0,0:n,:])*dts*2 # rain equation, no advection
      #r[2,1:n+1,:] = r[0,1:n+1,:] - alpha*dts*2.0*r[1,1:n+1,:]-2.0*beta_new[1:n+1,:]*(dts/dx)*(u[1,2:n+2,:]-u[1,1:n+1,:]) + (kr/(dx*dx))*(r[0,2:n+2,:] - 2.0*r[0,1:n+1,:] + r[0,0:n,:])*dts*2  - (dts/(dx))*(u[1,1:n+1,:])*(r[1,2:n+2,:]-r[1,0:n,:])  # other advection term 
     
      u[2,0,:]   = u[2,n,:]                                                  # boundary conditions
      u[2,n+1,:] = u[2,1,:]
      h[2,0,:]   = h[2,n,:]
      h[2,n+1,:] = h[2,1,:]
      r[2,0,:]   = r[2,n,:]
      r[2,n+1,:] = r[2,1,:]
      r[2,:,:] = np.where(r[2,:,:]<r_threshold,0,r[2,:,:])

      d = filter_parameter*.5*(u[2,:,:] - 2.*u[1,:,:] + u[0,:,:])            # RAW filter. Accounts for the growing computational mode.
      u[0,:,:] = u[1,:,:] + alpha_filt*d
      u[1,:,:] = u[2,:,:] - (1-alpha_filt)*d                     
      d = filter_parameter*.5*(h[2,:,:] - 2.*h[1,:,:] + h[0,:,:])
      h[0,:,:] = h[1,:,:] + alpha_filt*d
      h[1,:,:] = h[2,:,:] - (1-alpha_filt)*d
      d = filter_parameter*.5*(r[2,:,:] - 2.*r[1,:,:] + r[0,:,:])
      r[0,:,:] = r[1,:,:] + alpha_filt*d
      r[1,:,:] = r[2,:,:] - (1-alpha_filt)*d

      d = filter_parameter*.5*(phi_c[:,2,:] - 2.*phi_c[:,1,:] + phi_c[:,0,:])
      phi_c[:,0,:] = phi_c[:,1,:] + alpha_filt*d
      phi_c[:,1,:] = phi_c[:,2,:] + (1-alpha_filt)*d

    y[0:n,:] = u[2,1:n+1,:]
    y[n:2*n,:] = h[2,1:n+1,:]
    y[2*n:3*n,:] = r[2,1:n+1,:] 
    y_phi_c[:,:,:] = phi_c[:,[1,2],:] 
    y_phi[0:n,:] = phi[1:n+1,:]

    return y, y_phi_c, y_phi

def sweqssamp__(x,phi_data,phi_c_data,nsteps,time_amp,n,n_ens,mindex,nindex,time,batch_number,members_per_batch): # ss for selected saving 
    """
    x is the concatenated state of size 3*n. This is the relevant model input and ouput.
    nsteps is the number of timesteps
    n is the total grid size
    n_ens is the number of ensemble members
    mindex = mountain index - 1: bell shaped mountain
    nindex = noise index - 1: stochastic Gaussian noise throughout grid domain
    time indicates whether on first time step. Used to install orography and keep count for random seed
    amp in name indicates time dependence on amplitude of S_rad and noise
    """  
    phi = np.zeros((n+2,n_ens))                                 # geopotential
    phi_c = np.zeros((n,3,n_ens))                                 # constant geopotential to use when h>h_c
    beta_new = np.zeros((n+2,n_ens))                            # returns matrix of n+2 columns. Shape (n+2)
    harray = np.zeros((n+2,n_ens))                              # array which will be able to alter the depth of fluid
    y = np.zeros((3*n,n_ens))                                   # final outputting array
    y_phi_c = np.zeros((n,2,n_ens))                               # final outputting arrays
    y_phi = np.zeros((n,n_ens))
    array = np.zeros((n,n_ens))                                 # to be used in calculation
    x_for = np.arange(0,n,1)
    #no = []
    #ra = []

    ku = 3400                                                   # diffusion constant for u, usually 2000
    kh = 1400                                                   # diffusion coefficient for depth h, usually 6000
    
    if mindex == 1:                                             # bell shaped mountain ridge
      ku = 3000                                                 # diffusion constant for u
      kh = ku                                                   # diffusion coefficient for depth h
      for i in range(n):
        harray[i+1,:] = amp*hw*hw/((float(i)-mu)*(float(i)-mu)+hw*hw)   # every ensemble member has same mountain

    if mindex == 2:                                             # k**(-1) power spectrum
      x_sin = np.arange(0,n)                                    # mindex #1 and #2 need work and checking
      s = 0.2+np.sin((x_sin-100)/500)
      np.random.seed(2)
      data = np.random.rand(n)-0.5                              # data centred around zero
      ps = np.abs(np.fft.fft(data,norm='ortho'))**2             # calculating the power spectrum
      for i in range(n_ens):
        harray[1:n+1,i]=ps*s

    harray[0,:] = harray[n,:]
    harray[n+1,:] = harray[1,:]

    u = np.zeros((3,n+2,n_ens))                                             # horizontal velocity returns matrix of 3 row and n+2 columns 
    h = np.zeros((3,n+2,n_ens))                                             # depth, returns matrix of 3 row and n+2 columns 
    r = np.zeros((3,n+2,n_ens))                                             # rain, returns matrix of 3 row and n+2 columns 
    u[0,1:n+1,:],u[1,1:n+1,:] = x[0:n,:],x[0:n,:]                           # filling in the middle of the u,h and r matrixes with input 
    h[0,1:n+1,:],h[1,1:n+1,:] = x[n:2*n,:],x[n:2*n,:]
    r[0,1:n+1,:],r[1,1:n+1,:] = x[2*n:3*n,:],x[2*n:3*n,:]
    
    u[1,0,:],u[0,0,:] = u[0,n,:],u[0,n,:]                                   # filling in the 2 outer columns of the u, h and r matrixes
    u[1,n+1,:],u[0,n+1,:] = u[0,1,:],u[0,1,:]                               # for boundary conditions from which the model effectively works. 
    h[1,0,:],h[0,0,:] = h[0,n,:],h[0,n,:]
    h[1,n+1,:],h[0,n+1,:] =  h[0,1,:],h[0,1,:]
    r[1,0,:],r[0,0,:] = r[0,n,:],r[0,n,:]
    r[1,n+1,:],r[0,n+1,:] = r[0,1,:],r[0,1,:]

    phi[1:n+1,:] = phi_data
    phi_c[:,:2,:] = phi_c_data
    
    phi[0,:] = phi[n,:] # do periodicness as with u,h and r above
    phi[n+1,:] = phi[1,:]

    if time == 0:                                                           # intall ridge - only happens at first time step
      h[:,:,:] = h[:,:,:] - harray[:,:]     

    for it in range(nsteps):                                                # loop over the given model timestep. Saves output after nsteps
      if nindex == 1:
        noise_seed = time+it
        time_amp_noise = time_amp+it                                   # counts the timestep for the noise
        u, noiseamp_output = noiseamp(u,n,n_ens,n_array,noise_seed,batch_number,members_per_batch,time_amp_noise)                             # returns the u matrix with the 2nd row with stochastic noise across the ensemble
        u[1,0,:] = u[1,n,:]
        u[1,n+1,:] = u[1,1,:]
       # no = np.append(no,noiseamp_output)

      array = np.where(h[1,1:n+1,:]<38.4,0,(u[1,2:n+2,:] - u[1,1:n+1,:])/dx)
      array = np.where(array[:,:]>0,0,array[:,:])
      con_av = np.mean(array[:,:],axis=0)
      con_av = np.nan_to_num(con_av)

      S_conv = gamma_2 * beta * con_av[:]

      rad_amp = 0.00005*np.cos((2*freq*np.pi/(tf*t_whole))*(time_amp_noise)) 
      rad_amp = np.where(rad_amp>0,0,rad_amp)

      S_for_slowlyvarying = for_amp*np.sin(((2*freq*np.pi)/n)*(x_for)-np.int32(time_amp_noise/time_for_gp))+for_loc
      S_for_slowlyvarying = np.where(S_for_slowlyvarying<0,0,S_for_slowlyvarying)
      S_for_slowlyvarying = np.broadcast_to(S_for_slowlyvarying,(n_ens,n))
      S_for_slowlyvarying = S_for_slowlyvarying.T

      phi_c[:,2,:] = phi_c[:,0,:] + 2*dts*(rad_amp - S_for_slowlyvarying + S_conv)    # S_for_2 relates to free-run. Clarify how to do DA. 

      if mindex == 1 or mindex==2:
        phi[1:n+1,:] = np.where( h[1,1:n+1,:]+harray[1:n+1,:] > h_cloud , phi_c[:,2,:]+g*harray[1:n+1,:], g*h[1,1:n+1,:] ) # if condition met, return phic in h matrix. If not, return g*h thing which would be the normal geopotential below the first threshold
      else:
        phi[1:n+1,:] = np.where( h[1,1:n+1,:] > h_cloud , phi_c[:,2,:], g*h[1,1:n+1,:] )

      phi[0,:]   = phi[n,:]                                             # boundary conditions for geopotential
      phi[n+1,:] = phi[1,:]
      phi[:,:] = phi[:,:] + gamma * r[1,:,:] 

      # shallow water equations =D
      u[2,1:n+1,:] = u[0,1:n+1,:] - (dts/(2*dx))*(u[1,2:n+2,:]**2 - u[1,0:n,:]**2) - (2*dts/dx)*(phi[1:n+1,:]-phi[0:n,:])  + (ku/(dx*dx))*(u[0,2:n+2,:] - 2*u[0,1:n+1,:] + u[0,0:n,:])*dts*2     # momentum equation  # fixed diffusion term
      h[2,1:n+1,:] = h[0,1:n+1,:] - (dts/dx)*(u[1,2:n+2,:]*(h[1,1:n+1,:]+h[1,2:n+2,:]) - u[1,1:n+1,:]*(h[1,0:n,:]+h[1,1:n+1,:])) + (kh/(dx*dx))*(h[0,2:n+2,:] - 2*h[0,1:n+1,:] + h[0,0:n,:])*dts*2   # continuity equation  # fixed diffusion term

      if mindex == 1 or mindex==2:
        mask = np.logical_and(h[1,1:n+1,:]+harray[1:n+1,:] > h_rain, u[1,2:n+2,:]-u[1,1:n+1,:] < 0)    # conditions for rain
      else:
        mask = np.logical_and(h[1,1:n+1,:] > h_rain, u[1,2:n+2,:]-u[1,1:n+1,:] < 0)   

      beta_new[1:n+1,:] = np.where( mask, beta , 0 )
      r[2,1:n+1,:] = r[0,1:n+1,:] - (dts/(2*dx))*(u[1,2:n+2,:]+u[1,1:n+1,:])*(r[1,2:n+2,:]-r[1,0:n,:]) - alpha*dts*2.0*r[1,1:n+1,:]-2.0*beta_new[1:n+1,:]*(dts/dx)*(u[1,2:n+2,:]-u[1,1:n+1,:]) + (kr/(dx*dx))*(r[0,2:n+2,:] - 2.0*r[0,1:n+1,:] + r[0,0:n,:])*dts*2 # rain equation  # with advection # always peaky
      #r[2,1:n+1,:] = r[0,1:n+1,:] - alpha*dts*2.0*r[1,1:n+1,:]-2.0*beta_new[1:n+1,:]*(dts/dx)*(u[1,2:n+2,:]-u[1,1:n+1,:]) + (kr/(dx*dx))*(r[0,2:n+2,:] - 2.0*r[0,1:n+1,:] + r[0,0:n,:])*dts*2 # rain equation, no advection
      #r[2,1:n+1,:] = r[0,1:n+1,:] - alpha*dts*2.0*r[1,1:n+1,:]-2.0*beta_new[1:n+1,:]*(dts/dx)*(u[1,2:n+2,:]-u[1,1:n+1,:]) + (kr/(dx*dx))*(r[0,2:n+2,:] - 2.0*r[0,1:n+1,:] + r[0,0:n,:])*dts*2  - (dts/(dx))*(u[1,1:n+1,:])*(r[1,2:n+2,:]-r[1,0:n,:])  # other advection term 
     
      u[2,0,:]   = u[2,n,:]                                                  # boundary conditions
      u[2,n+1,:] = u[2,1,:]
      h[2,0,:]   = h[2,n,:]
      h[2,n+1,:] = h[2,1,:]
      r[2,0,:]   = r[2,n,:]
      r[2,n+1,:] = r[2,1,:]
      r[2,:,:] = np.where(r[2,:,:]<r_threshold,0,r[2,:,:])

      d = filter_parameter*.5*(u[2,:,:] - 2.*u[1,:,:] + u[0,:,:])            # RAW filter. Accounts for the growing computational mode.
      u[0,:,:] = u[1,:,:] + alpha_filt*d
      u[1,:,:] = u[2,:,:] - (1-alpha_filt)*d                     
      d = filter_parameter*.5*(h[2,:,:] - 2.*h[1,:,:] + h[0,:,:])
      h[0,:,:] = h[1,:,:] + alpha_filt*d
      h[1,:,:] = h[2,:,:] - (1-alpha_filt)*d
      d = filter_parameter*.5*(r[2,:,:] - 2.*r[1,:,:] + r[0,:,:])
      r[0,:,:] = r[1,:,:] + alpha_filt*d
      r[1,:,:] = r[2,:,:] - (1-alpha_filt)*d

      d = filter_parameter*.5*(phi_c[:,2,:] - 2.*phi_c[:,1,:] + phi_c[:,0,:])
      phi_c[:,0,:] = phi_c[:,1,:] + alpha_filt*d
      phi_c[:,1,:] = phi_c[:,2,:] + (1-alpha_filt)*d

    y[0:n,:] = u[2,1:n+1,:]
    y[n:2*n,:] = h[2,1:n+1,:]
    y[2*n:3*n,:] = r[2,1:n+1,:] 
    y_phi_c[:,:] = phi_c[:,[1,2],:] 
    y_phi[0:n,:] = phi[1:n+1,:]
    print('time_amp_noise: '+str(time_amp_noise))

    #np.savetxt(str(save_directory)+'no_'+str(time_amp),no)
    #np.savetxt(str(save_directory)+'ra'+str(time_amp),ra)

    return y, y_phi_c, y_phi


def Initialize__(n,nens,mindex,nindex):
  '''
  Initialise x array input for sweq function when initialising DA. 
  n - grid size 
  nens - number of ensembles. Note!: no multiplying factor here since 'Initialize' used for DA only.
  '''
  state = np.zeros((3*n,nens))
  state[0:n] = np.zeros((n,nens))+0.                          # velocity 0m/s
  state[n:2*n] = np.zeros((n,nens))+38.                       # fluid depth 38m
  phi_data = np.zeros((n,nens))
  phi_c_data = np.zeros((n,2,nens))
  phi_data = phi_data + 380
  phi_c_data = phi_c_data + 380 # 899.96

  state, y_phi_c, y_phi = sweq__(state,phi_data,phi_c_data,1000,n,nens,mindex,nindex,time=0,batch_number=0,members_per_batch=0)      # state to begin with. 1000 nsteps taken
  return state,  y_phi_c, y_phi 

def Initializenew__(n,nens,mindex,nindex):
  '''
  Initialise x array input for sweq function when initialising DA. 
  n - grid size 
  nens - number of ensembles. Note!: no multiplying factor here since 'Initialize' used for DA only.
  new - no simulation carried out. 
  '''
  state = np.zeros((3*n,nens))
  state[0:n] = np.zeros((n,nens))+0.                          # velocity 0m/s
  state[n:2*n] = np.zeros((n,nens))+38.                       # fluid depth 38m
  y_phi = np.zeros((n,nens)) + 380
  y_phi_c = np.zeros((n,2,nens)) + 380

  #state, y_phi_c, y_phi = sweq__(state,phi_data,phi_c_data,1000,n,nens,mindex,nindex,time=0,batch_number=0,members_per_batch=0)      # state to begin with. 1000 nsteps taken
  return state,  y_phi_c, y_phi # dimensions different to Initialise function - no time component

def noise(u,n,n_ens,n_array,noise_seed,batch_number,members_per_batch):   
    ''' 
    Random stochastic noise added to the velocity. Mimics turbulence in the boundary layer. 
    n_array tells model where noise is located. narray.shape = (number of perturbation areas, features of pert area)
    features of pert area = [half width, amplitude, start, end of noise field region]
    iteration added so 'random' noise is seeded but different on each iteration cycle
    '''
    unoise = np.zeros((2*n,n_ens))
    mu = float((n+1.0)/2.0)                                           # center of noise

    for area in range(len(n_array[:,0])):                             # loop goes over all areas of perturbations 
      sig = n_array[area,0]                                           # sigma, half width of noise field
      amp = n_array[area,1]                                           # amplitude of noise field
      d = np.array(range(n+1))
      z = (1.0/(sig*np.sqrt(2.0*np.pi)))*np.exp(-0.5*((d-mu)/sig)**2)
      zsum = z[1:n+1]-z[0:n]                                          # zsum will look like convergence
      zsum = amp*zsum/max(zsum)                                       # normalising

      for e in range(n_ens):                                          # adds this noise to all ensembles
        start = n_array[area,2]
        end = n_array[area,3]
        e_update = e + (batch_number*members_per_batch)               # So each batch doesn't recieve the same forcing
        noise_seed_chosen = ((noise_seed+e_update)*(e_update+2))-(noise_seed+e_update*3)   # choosing the random seed so that each ensemble member is perturbed differently, and at different times differently also. But same seed chosen for all orographic features. 
        while noise_seed_chosen>=2**32:
            noise_seed_chosen = noise_seed_chosen-(2**32)             # returns error if seed above this value
        np.random.seed(noise_seed_chosen)  
        pos = np.random.randint(start,end)                            # each ensemble will have own random noise added within bounds created
        unoise[pos:pos+n,e] = unoise[pos:pos+n,e] + zsum
        unoise[0:n,e] = unoise[0:n,e] + unoise[n:n+n,e]
        unoise[n:n+n,e] = 0 
        u[1,1:n+1,e] = u[1,1:n+1,e] + unoise[0:n,e]
    return u

def noiseamp(u,n,n_ens,n_array,noise_seed,batch_number,members_per_batch,time_amp_noise):   
    ''' 
    Random stochastic noise added to the velocity. Mimics turbulence in the boundary layer. 
    n_array tells model where noise is located. narray.shape = (number of perturbation areas, features of pert area)
    features of pert area = [half width, amplitude, start, end of noise field region]
    iteration added so 'random' noise is seeded but different on each iteration cycle
    amp refers to time dependence on the amplitude of the noise
    '''
    unoise = np.zeros((2*n,n_ens))
    mu = float((n+1.0)/2.0)                                           # center of noise
    amp_cos = -np.cos((2*freq*np.pi/(tf*t_whole))*time_amp_noise)        # time dependent amplitude
    amp_cos = np.where(amp_cos<0,0,amp_cos)
    #amp_cos = np.where(amp_cos<0,0.5,amp_cos+0.5)              # use if want to add in (unrealistic) base noise

    for area in range(len(n_array[:,0])):                             # loop goes over all areas of perturbations 
      sig = n_arrayamp[area,0]                                           # sigma, half width of noise field
      amp = n_arrayamp[area,1]                                           # amplitude of noise field
      d = np.array(range(n+1))
      z = (1.0/(sig*np.sqrt(2.0*np.pi)))*np.exp(-0.5*((d-mu)/sig)**2)
      zsum = z[1:n+1]-z[0:n]                                          # zsum will look like convergence
      zsum = amp_cos*amp*zsum/max(zsum)                                       # normalising

      for e in range(n_ens):                                          # adds this noise to all ensembles
        start = n_arrayamp[area,2]
        end = n_arrayamp[area,3]
        e_update = e + (batch_number*members_per_batch)               # So each batch doesn't recieve the same forcing
        noise_seed_chosen = ((noise_seed+e_update)*(e_update+2))-(noise_seed+e_update*3)   # choosing the random seed so that each ensemble member is perturbed differently, and at different times differently also. But same seed chosen for all orographic features. 
        while noise_seed_chosen>=2**32:
            noise_seed_chosen = noise_seed_chosen-(2**32)             # returns error if seed above this value
        np.random.seed(noise_seed_chosen)  
        pos = np.random.randint(start,end)                            # each ensemble will have own random noise added within bounds created
        unoise[pos:pos+n,e] = unoise[pos:pos+n,e] + zsum
        unoise[0:n,e] = unoise[0:n,e] + unoise[n:n+n,e]
        unoise[n:n+n,e] = 0 
        u[1,1:n+1,e] = u[1,1:n+1,e] + unoise[0:n,e]
    
    return u, amp_cos

def model_initss__(start_time,batch_number,members_per_batch):   # ss for selected saving 
  '''
  Running the SWE with the IC from DA_init_vis_dt__)
  Tracks time
  Multipy says by what factor to increase the ensemble size
  Creates 'actual_truth' run which is the truth from the DA initialisation continued. So it really has a total ensemble size of (nens*multiply)+1
  If multiply = True, then enemble multiplied. Don't want this if continuing run. 
  Need to change actual_truth and DA_time dependent on which DA used. 
  '''
  start_time_model = datetime.now()
  if start_time==0:
      SWMDA = Dataset(str(save_directory)+'SWMDA.nc','r')                # read DA netCDF
      SWM = Dataset(str(save_directory)+'SWM0_'+str(batch_number)+'_'+str(members_per_batch)+'.nc','w')    # want to create a new file, don't use same file as DA anymore!!!
  else:
      OLD = Dataset(str(save_directory)+'SWM'+str(start_time-t)+'_'+str(batch_number)+'_'+str(members_per_batch)+'.nc','r')
      SWM = Dataset(str(save_directory)+'SWM'+str(start_time)+'_'+str(batch_number)+'_'+str(members_per_batch)+'.nc','w',format='NETCDF4') 
      
  free_run = SWM.createGroup('free_run')                                 # create free_run group
  no_ens = free_run.createDimension('number of members',members_per_batch) # create dimensions
  time = free_run.createDimension('time',None)
  total_data = free_run.createDimension('total data',3*n)
  one = free_run.createDimension('one',1)
  two = free_run.createDimension('twoS',2)
  grid_points = free_run.createDimension('grid_pointS',n)

  truth = free_run.createVariable('truth','f8',('total data','number of members','time')) # create variables
  truth_phi_c = free_run.createVariable('truth_phi_c','f8',('grid_pointS','twoS','number of members','time'))
  truth_phi = free_run.createVariable('truth_phi','f8',('grid_pointS','number of members','time'))
  actual_truth = free_run.createVariable('actual_truth','f8',('total data','one','time'))
  actual_truth_phi_c = free_run.createVariable('actual_truth_phi_c','f8',('grid_pointS','twoS','one','time'))
  actual_truth_phi = free_run.createVariable('actual_truth_phi','f8',('grid_pointS','one','time'))

  DA_time = (1+ncyc)*(dte)+1000         
  if start_time==0:             
      actual_truth[:,:,0] = ma.getdata(SWMDA['/DA/truth_plot'][:,:,(ncyc*dte)-1])        # picking up the truth from the DA initialisation
      actual_truth_phi_c[:,:,:,0] = ma.getdata(SWMDA['/DA/truth_plot_phi_c'][:,:,:,(ncyc*dte)-1])
      actual_truth_phi[:,:,0] = ma.getdata(SWMDA['/DA/truth_plot_phi'][:,:,(ncyc*dte)-1])
  else:
      actual_truth[:,:,0] = ma.getdata(OLD['/free_run/actual_truth'][:,:,t-1])
      actual_truth_phi_c[:,:,:,0] = ma.getdata(OLD['/free_run/actual_truth_phi_c'][:,:,:,t-1])
      actual_truth_phi[:,:,0] = ma.getdata(OLD['/free_run/actual_truth_phi'][:,:,t-1])

  for i in range(t):
    print('Computing the '+str(i)+' of t now.')
    if i == 0:
      actual_truth[:,:,0], actual_truth_phi_c[:,:,:,0], actual_truth_phi[:,:,0]= sweqss__(actual_truth[:,:,0],actual_truth_phi[:,:,0],actual_truth_phi_c[:,:,:,0],nsteps=tf,n=n,n_ens=1,mindex=mindex,nindex=nindex,time=(start_time*tf)+DA_time,batch_number=0,members_per_batch=0) 
    else:
      actual_truth[:,:,i],actual_truth_phi_c[:,:,:,i],actual_truth_phi[:,:,i] = sweqss__(actual_truth[:,:,i-1],actual_truth_phi[:,:,i-1],actual_truth_phi_c[:,:,:,i-1],nsteps=tf,n=n,n_ens=1,mindex=mindex,nindex=nindex,time=(start_time*tf)+DA_time+(i*tf),batch_number=0,members_per_batch=0) # integrating this 'real' truth further in time for model simulation

  if start_time ==0:                                                                                                          # first time step of truth given by DA initialisation
      truth[:,:,0] = np.repeat(ma.getdata(SWMDA['/DA/analysis'][ncyc-1,:,batch_number*split:(batch_number+1)*split]),repeats=multiply,axis=1)          # change this line depending on how splitting up analysis!                 # method depends on whether or not you wish to expand the ensemble size
      truth_phi_c[:,:,:,0] = np.repeat(ma.getdata(SWMDA['/DA/background_phi_c'][:,:,batch_number*split:(batch_number+1)*split,(ncyc*dte)-1]),repeats=multiply,axis=1) 
      truth_phi[:,:,0] = np.repeat(ma.getdata(SWMDA['/DA/background_phi'][:,batch_number*split:(batch_number+1)*split,(ncyc*dte)-1]),repeats=multiply,axis=1) 
      if len(truth[1,:,1]) == members_per_batch:
        print('split correctly')
      else:
        print('not split correctly!!!!')
  else:
      truth[:,:,0] = ma.getdata(OLD['/free_run/truth'][:,:,t-1])
      truth_phi_c[:,:,:,0] = ma.getdata(OLD['/free_run/truth_phi_c'][:,:,:,t-1])
      truth_phi[:,:,0] = ma.getdata(OLD['/free_run/truth_phi'][:,:,t-1])

  for i in range(t):
    if i == 0:
      truth[:,:,0], truth_phi_c[:,:,:,0],truth_phi[:,:,0] = sweqss__(truth[:,:,0], truth_phi[:,:,0],truth_phi_c[:,:,:,0],nsteps=tf,n=n,n_ens=members_per_batch,mindex=mindex,nindex=nindex,time=(start_time*tf)+DA_time,batch_number=batch_number,members_per_batch=members_per_batch) 
    else:
      truth[:,:,i], truth_phi_c[:,:,:,i],truth_phi[:,:,i] = sweqss__(truth[:,:,i-1], truth_phi[:,:,i-1], truth_phi_c[:,:,:,i-1],nsteps=tf,n=n,n_ens=members_per_batch,mindex=mindex,nindex=nindex,time=(start_time*tf)+DA_time+(i*tf),batch_number=batch_number,members_per_batch=members_per_batch)     # generate model simulation, store in 'truth'

  time = datetime.now()-start_time_model
  print('time for model run with nens= ',members_per_batch,' equals= ', time)
  SWMDA.close() 
  SWM.close() 
  return 

def model_initssamp__(start_time,batch_number,members_per_batch):   # ss for selected saving 
  '''
  Running the SWE with the IC from DA_init_vis_dt__)
  Tracks time
  Multipy says by what factor to increase the ensemble size
  Creates 'actual_truth' run which is the truth from the DA initialisation continued. So it really has a total ensemble size of (nens*multiply)+1
  If multiply = True, then enemble multiplied. Don't want this if continuing run. 
  Now using sweq function with amp. 
  Edited for new DA amp initialisation function
  '''
  start_time_model = datetime.now()
  if start_time==0:
      SWMDA = Dataset(str(save_directory)+'SWMDA.nc','r')                # read DA netCDF
      SWM = Dataset(str(save_directory)+'SWM0_'+str(batch_number)+'_'+str(members_per_batch)+'.nc','w')    # want to create a new file, don't use same file as DA anymore!!!
  else:
      OLD = Dataset(str(save_directory)+'SWM'+str(start_time-t)+'_'+str(batch_number)+'_'+str(members_per_batch)+'.nc','r')
      SWM = Dataset(str(save_directory)+'SWM'+str(start_time)+'_'+str(batch_number)+'_'+str(members_per_batch)+'.nc','w',format='NETCDF4') 
      
  free_run = SWM.createGroup('free_run')                                 # create free_run group
  no_ens = free_run.createDimension('number of members',members_per_batch) # create dimensions
  time = free_run.createDimension('time',None)
  total_data = free_run.createDimension('total data',3*n)
  one = free_run.createDimension('one',1)
  two = free_run.createDimension('twoS',2)
  grid_points = free_run.createDimension('grid_pointS',n)

  truth = free_run.createVariable('truth','f8',('total data','number of members','time')) # create variables
  truth_phi_c = free_run.createVariable('truth_phi_c','f8',('grid_pointS','twoS','number of members','time'))
  truth_phi = free_run.createVariable('truth_phi','f8',('grid_pointS','number of members','time'))
  actual_truth = free_run.createVariable('actual_truth','f8',('total data','one','time'))
  actual_truth_phi_c = free_run.createVariable('actual_truth_phi_c','f8',('grid_pointS','twoS','one','time'))
  actual_truth_phi = free_run.createVariable('actual_truth_phi','f8',('grid_pointS','one','time'))

  DA_time = ncyc*dte       
  if start_time==0:             
      actual_truth[:,:,0] = ma.getdata(SWMDA['/DA/truth_plot'][:,:,ncyc-1])        # picking up the truth from the DA initialisation
      actual_truth_phi_c[:,:,:,0] = ma.getdata(SWMDA['/DA/truth_plot_phi_c'][:,:,:,ncyc-1])
      actual_truth_phi[:,:,0] = ma.getdata(SWMDA['/DA/truth_plot_phi'][:,:,ncyc-1])
  else:
      actual_truth[:,:,0] = ma.getdata(OLD['/free_run/actual_truth'][:,:,t-1])
      actual_truth_phi_c[:,:,:,0] = ma.getdata(OLD['/free_run/actual_truth_phi_c'][:,:,:,t-1])
      actual_truth_phi[:,:,0] = ma.getdata(OLD['/free_run/actual_truth_phi'][:,:,t-1])

  for i in range(t):
    print('Computing the '+str(i)+' of t now.')
    if i == 0:
      actual_truth[:,:,0], actual_truth_phi_c[:,:,:,0], actual_truth_phi[:,:,0]= sweqssamp__(actual_truth[:,:,0],actual_truth_phi[:,:,0],actual_truth_phi_c[:,:,:,0],nsteps=tf,time_amp=start_time*tf,n=n,n_ens=1,mindex=mindex,nindex=nindex,time=(start_time*tf)+DA_time,batch_number=0,members_per_batch=0) 
    else:
      actual_truth[:,:,i],actual_truth_phi_c[:,:,:,i],actual_truth_phi[:,:,i] = sweqssamp__(actual_truth[:,:,i-1],actual_truth_phi[:,:,i-1],actual_truth_phi_c[:,:,:,i-1],nsteps=tf,time_amp=(i+start_time)*tf,n=n,n_ens=1,mindex=mindex,nindex=nindex,time=(start_time*tf)+DA_time+(i*tf),batch_number=0,members_per_batch=0) # integrating this 'real' truth further in time for model simulation

  if start_time ==0:                                                                                                          # first time step of truth given by DA initialisation
      truth[:,:,0] = np.repeat(ma.getdata(SWMDA['/DA/analysis'][ncyc-1,:,batch_number*split:(batch_number+1)*split]),repeats=multiply,axis=1)          # change this line depending on how splitting up analysis!                 # method depends on whether or not you wish to expand the ensemble size
      truth_phi_c[:,:,:,0] = np.repeat(ma.getdata(SWMDA['/DA/background_phi_c'][:,:,batch_number*split:(batch_number+1)*split,ncyc-1]),repeats=multiply,axis=1) 
      truth_phi[:,:,0] = np.repeat(ma.getdata(SWMDA['/DA/background_phi'][:,batch_number*split:(batch_number+1)*split,ncyc-1]),repeats=multiply,axis=1) 
      if len(truth[1,:,1]) == members_per_batch:
        print('split correctly')
      else:
        print('not split correctly!!!!')
  else:
      truth[:,:,0] = ma.getdata(OLD['/free_run/truth'][:,:,t-1])
      truth_phi_c[:,:,:,0] = ma.getdata(OLD['/free_run/truth_phi_c'][:,:,:,t-1])
      truth_phi[:,:,0] = ma.getdata(OLD['/free_run/truth_phi'][:,:,t-1])

  for i in range(t):
    if i == 0:
      truth[:,:,0], truth_phi_c[:,:,:,0],truth_phi[:,:,0] = sweqssamp__(truth[:,:,0], truth_phi[:,:,0],truth_phi_c[:,:,:,0],nsteps=tf,time_amp=start_time*tf,n=n,n_ens=members_per_batch,mindex=mindex,nindex=nindex,time=(start_time*tf)+DA_time,batch_number=batch_number,members_per_batch=members_per_batch) 
    else:
      truth[:,:,i], truth_phi_c[:,:,:,i],truth_phi[:,:,i] = sweqssamp__(truth[:,:,i-1], truth_phi[:,:,i-1], truth_phi_c[:,:,:,i-1],nsteps=tf,time_amp=(i+start_time)*tf,n=n,n_ens=members_per_batch,mindex=mindex,nindex=nindex,time=(start_time*tf)+DA_time+(i*tf),batch_number=batch_number,members_per_batch=members_per_batch)     # generate model simulation, store in 'truth'

  time = datetime.now()-start_time_model
  print('time for model run with nens= ',members_per_batch,' equals= ', time)
  SWMDA.close() 
  SWM.close() 
  return 

# def DA_init_vis_dt__(DA_method):
#   '''
#   Initialising the model with DA using the DA_method(usually EnKF). Goes through ncyc cycles before outoutting the array model_start which are the IC. 
#   Also returns rmse error of the background and analysis, bk_error and an_error
#   Tracks time for initialisation
#   This function allows to define a different dt compared to the model run. dte is defined in the constants.
#   '''
#   SWM = Dataset(str(save_directory)+'SWMDA.nc','w',format='NETCDF4')                # set up netCDF file for saving output
#   SWM.description = 'shallow water model'
#   SWM.history = 'Created ' + str(datetime.now())

#   one = SWM.createDimension('oneS',1)
#   two = SWM.createDimension('two',2)
#   four = SWM.createDimension('four',4)

#   grid_point = SWM.createDimension('grid_points',1)           # inputing contants to nc file
#   grid_points = SWM.createVariable('n','i8',('grid_points',))
#   grid_points[:] = n
#   dte_const = SWM.createDimension('model time step for DA initialisation',1)
#   dte_consts = SWM.createVariable('dte','i8',('model time step for DA initialisation',))
#   dte_consts[:] = dte
#   t_const = SWM.createDimension('length of the model simulation',1)
#   t_consts = SWM.createVariable('t','i8',('length of the model simulation',))
#   t_consts[:] = t
#   multiply_const = SWM.createDimension('factor by which to increase ensemble size after DA initialisation',1)
#   multiply_consts = SWM.createVariable('multiply','i8',('factor by which to increase ensemble size after DA initialisation',))
#   multiply_consts[:] = multiply
#   sigma_const = SWM.createDimension('standard deviation of observation error',3)
#   sigma_consts = SWM.createVariable('sigma','f8',('standard deviation of observation error',))
#   sigma_consts[:] = sig
#   loc_const = SWM.createDimension('localisation radius',1)
#   loc_consts = SWM.createVariable('localisation radius','i8',('localisation radius',))
#   loc_consts[:] = loc_radius
#   mindex_const = SWM.createDimension('mindex',1)
#   mindex_consts = SWM.createVariable('mountain index','u8',('mindex',))
#   mindex_consts[:] = mindex
#   nindex_const = SWM.createDimension('nindex',1)
#   nindex_consts = SWM.createVariable('noise index','u8',('nindex',))
#   nindex_consts[:] = nindex
#   n_array_consts = SWM.createVariable('noise array','f8',('oneS','four',))
#   n_array_consts[:] = n_array
   
#   DA = SWM.createGroup('DA')                                  # create DA group
#   total_data = DA.createDimension('total data',3*n)           # create dimensions
#   no_ens = DA.createDimension('number of members',nens)
#   no_cyc = DA.createDimension('number of DA cycles',ncyc)
#   no_cyc1 = DA.createDimension('ncyc + 1',ncyc+1)
#   one = DA.createDimension('one',1)
#   truth_data_points = DA.createDimension('truth data points',dte*(ncyc+1))
#   init_length = DA.createDimension('init length',1000)
#   two = DA.createDimension('twoS',2)
#   grid_points = SWM.createDimension('grid_pointS',n)

#   analysis = DA.createVariable('analysis','f8',('number of DA cycles','total data','number of members'))  # create variables
#   background = DA.createVariable('background','f8',('total data','number of members','truth data points'))
#   background_phi_c = DA.createVariable('background_phi_c','f8',('grid_pointS','twoS','number of members','truth data points'))
#   background_phi = DA.createVariable('background_phi','f8',('grid_pointS','number of members','truth data points'))
#   observation = DA.createVariable('observation','f8',('number of DA cycles','total data'))
#   truth_plot = DA.createVariable('truth_plot','f8',('total data','one','truth data points'))
#   truth_plot_phi_c = DA.createVariable('truth_plot_phi_c','f8',('grid_pointS','twoS','one','truth data points'))
#   truth_plot_phi = DA.createVariable('truth_plot_phi','f8',('grid_pointS','one','truth data points'))
#   truth_plot_init = DA.createVariable('truth_plot_init','f8',('total data','one','init length'))
#   truth_plot_phi_c_init = DA.createVariable('truth_plot_phi_c_init','f8',('grid_pointS','twoS','one','init length'))
#   truth_plot_phi_init = DA.createVariable('truth_plot_phi_init','f8',('grid_pointS','one','init length'))
#   background_init = DA.createVariable('background_init','f8',('total data','number of members','init length'))
#   background_init_phi_c = DA.createVariable('background_init_phi_c','f8',('grid_pointS','twoS','number of members','init length'))
#   background_init_phi = DA.createVariable('background_init_phi','f8',('grid_pointS','number of members','init length'))

#   start_time = datetime.now()
#   bg_error = {}                                               # background error dictionary
#   an_error = {}                                               # analysis error dictionary
#   for var in ['u','h','r']:
#     bg_error[var] = np.zeros((ncyc))                   
#     an_error[var] = np.zeros((ncyc))
   
#   truth_plot_init[:,:,:], truth_plot_phi_c_init[:,:,:,:], truth_plot_phi_init[:,:,:] = Initialize__(n,1,mindex,nindex)
#   truth_plot[:,:,0] = truth_plot_init[:,:,999]                # use last time step
#   truth_plot_phi_c[:,:,:,0] = truth_plot_phi_c_init[:,:,:,999]
#   truth_plot_phi[:,:,0] = truth_plot_phi_init[:,:,999]
#   truth_plot[:,:,1:dte], truth_plot_phi_c[:,:,:,1:dte], truth_plot_phi[:,:,1:dte]  = sweq__(truth_plot[:,:,0],truth_plot_phi[:,:,0],truth_plot_phi_c[:,:,:,0],nsteps=dte-1,n=n,n_ens=1,mindex=mindex,nindex=nindex,time=1000,batch_number=0,members_per_batch=0)                # generate the truth model simulation. Time parameter = 0 to install mountain since first iteration
#   truth = truth_plot[:,:,dte-1]
#   truth_phi_c = truth_plot_phi_c[:,:,:,dte-1]
#   truth_phi = truth_plot_phi[:,:,dte-1]
  
#   background_init[:,:,:], background_init_phi_c[:,:,:,:], background_init_phi[:,:,:] = Initialize__(n,nens,mindex,nindex)
#   background[:,:,0] = background_init[:,:,999]                # generate the initial ensemble, use last time step
#   background_phi_c[:,:,:,0] = background_init_phi_c[:,:,:,999]
#   background_phi[:,:,0] = background_init_phi[:,:,999]
    
#   background[:,:,1:dte], background_phi_c[:,:,:,1:dte], background_phi[:,:,1:dte] = sweq__(background[:,:,0],background_phi[:,:,0],background_phi_c[:,:,:,0],nsteps=dte-1,n=n,n_ens=nens,mindex=mindex,nindex=nindex,time=1000,batch_number=0,members_per_batch=0)             # propagate this ensemble. Time parameter = 0 to install mountain since forst iteration
#   model = background[:,:,dte-1]
#   model_phi_c = background_phi_c[:,:,:,dte-1]
#   model_phi = background_phi[:,:,dte-1]

#   for i in range(ncyc):                                       # go through ncyc DA cycles
#       print('now doing cycle: '+str(i))
#       for j,var in enumerate(['u','h','r']):                  # append background error 
#          bg_error[var][i] = rmse(truth[j*n:n*(j+1),0],np.mean(model[j*n:n*(j+1),:],axis=1))
#       truth_obs = truth.reshape(-1)                           # reshape makes shape of array (3*n) instead of (3*n,1). Needed for get_obs function
#       obs_ncyc = i 
#       obs = get_obs(truth=truth_obs,seed=obs_ncyc)            # create observations
#       observation[i,:] = obs
#       model = DA_method(obs,model,obs_ncyc)                   # compute analysis
#       analysis[i,:,:] = model
     
#       for j,var in enumerate(['u','h','r']):                  # append analysis error
#           an_error[var][i] = rmse(truth[j*n:n*(j+1),0],np.mean(model[j*n:n*(j+1),:],axis=1))

#       truth_plot[:,:,((i*dte)+dte):((i*dte)+(2*dte))], truth_plot_phi_c[:,:,:,((i*dte)+dte):((i*dte)+(2*dte))], truth_plot_phi[:,:,((i*dte)+dte):((i*dte)+(2*dte))] = sweq__(truth,truth_phi, truth_phi_c,nsteps=dte,n=n,n_ens=1,mindex=mindex,nindex=nindex,time=1000+((1+i)*dte),batch_number=0,members_per_batch=0)
#       truth = truth_plot[:,:,(2*dte)+(i*dte)-1]
#       truth_phi_c = truth_plot_phi_c[:,:,:,(2*dte)+(i*dte)-1]
#       truth_phi = truth_plot_phi[:,:,(2*dte)+(i*dte)-1]

#       background[:,:,((i*dte)+dte):((i*dte)+(2*dte))], background_phi_c[:,:,:,((i*dte)+dte):((i*dte)+(2*dte))], background_phi[:,:,((i*dte)+dte):((i*dte)+(2*dte))] = sweq__(model, model_phi, model_phi_c,nsteps=dte,n=n,n_ens=nens,mindex=mindex,nindex=nindex,time=1000+((1+i)*dte),batch_number=0,members_per_batch=0)
#       model = background[:,:,(2*dte)+(i*dte)-1]
#       model_phi_c = background_phi_c[:,:,:,(2*dte)+(i*dte)-1]
#       model_phi = background_phi[:,:,(2*dte)+(i*dte)-1]

#   pickle_obj = open(str(save_directory)+'dicts.an_error','wb')                  # saving an_error and bg_error in dictionary
#   pickle.dump(an_error,pickle_obj)
#   pickle_obj.close()
#   pickle_obj = open(str(save_directory)+'dicts.bg_error','wb')
#   pickle.dump(bg_error,pickle_obj)
#   pickle_obj.close()
      
#   time = datetime.now()-start_time
#   print('time for DA initialisation with nens= ',nens,' and ncyc= ',ncyc,' equals= ', time)

#   SWM.close()

#   return

def DA_init_vis_dtampss__(DA_method):
  '''
  Initialising the model with DA using the DA_method(usually EnKF). Goes through ncyc cycles before outoutting the array model_start which are the IC. 
  Also returns rmse error of the background and analysis, bk_error and an_error
  Tracks time for initialisation
  This function allows to define a different dt compared to the model run. dte is defined in the constants.
  This new updated (-amp__) version doesn't use an initial 1000 timestep initialisation and goes through the diurnal cycle. Also now use selected saving (s).
  '''
  SWM = Dataset(str(save_directory)+'SWMDA.nc','w',format='NETCDF4')                # set up netCDF file for saving output
  SWM.description = 'shallow water model'
  SWM.history = 'Created ' + str(datetime.now())

  one = SWM.createDimension('oneS',1)
  two = SWM.createDimension('two',2)
  four = SWM.createDimension('four',4)

  grid_point = SWM.createDimension('grid_points',1)           # inputing contants to nc file
  grid_points = SWM.createVariable('n','i8',('grid_points',))
  grid_points[:] = n
  dte_const = SWM.createDimension('model time step for DA initialisation',1)
  dte_consts = SWM.createVariable('dte','i8',('model time step for DA initialisation',))
  dte_consts[:] = dte
  t_const = SWM.createDimension('length of the model simulation',1)
  t_consts = SWM.createVariable('t','i8',('length of the model simulation',))
  t_consts[:] = t
  multiply_const = SWM.createDimension('factor by which to increase ensemble size after DA initialisation',1)
  multiply_consts = SWM.createVariable('multiply','i8',('factor by which to increase ensemble size after DA initialisation',))
  multiply_consts[:] = multiply
  sigma_const = SWM.createDimension('standard deviation of observation error',3)
  sigma_consts = SWM.createVariable('sigma','f8',('standard deviation of observation error',))
  sigma_consts[:] = sig
  loc_const = SWM.createDimension('localisation radius',1)
  loc_consts = SWM.createVariable('localisation radius','i8',('localisation radius',))
  loc_consts[:] = loc_radius
  mindex_const = SWM.createDimension('mindex',1)
  mindex_consts = SWM.createVariable('mountain index','u8',('mindex',))
  mindex_consts[:] = mindex
  nindex_const = SWM.createDimension('nindex',1)
  nindex_consts = SWM.createVariable('noise index','u8',('nindex',))
  nindex_consts[:] = nindex
  n_array_consts = SWM.createVariable('noise array','f8',('oneS','four',))
  n_array_consts[:] = n_array
   
  DA = SWM.createGroup('DA')                                  # create DA group
  total_data = DA.createDimension('total data',3*n)           # create dimensions
  no_ens = DA.createDimension('number of members',nens)
  no_cyc = DA.createDimension('number of DA cycles',ncyc)
  no_cyc1 = DA.createDimension('ncyc + 1',ncyc+1)
  one = DA.createDimension('one',1)
  two = DA.createDimension('twoS',2)
  grid_points = SWM.createDimension('grid_pointS',n)

  analysis = DA.createVariable('analysis','f8',('number of DA cycles','total data','number of members'))  # create variables
  background = DA.createVariable('background','f8',('total data','number of members','ncyc + 1'))  # background used to get analysis
  background_phi_c = DA.createVariable('background_phi_c','f8',('grid_pointS','twoS','number of members','ncyc + 1'))
  background_phi = DA.createVariable('background_phi','f8',('grid_pointS','number of members','ncyc + 1'))
  observation = DA.createVariable('observation','f8',('number of DA cycles','total data'))
  truth_plot = DA.createVariable('truth_plot','f8',('total data','one','ncyc + 1'))
  truth_plot_phi_c = DA.createVariable('truth_plot_phi_c','f8',('grid_pointS','twoS','one','ncyc + 1'))
  truth_plot_phi = DA.createVariable('truth_plot_phi','f8',('grid_pointS','one','ncyc + 1'))
  truth_plot_init = DA.createVariable('truth_plot_init','f8',('total data','one'))
  truth_plot_phi_c_init = DA.createVariable('truth_plot_phi_c_init','f8',('grid_pointS','twoS','one'))
  truth_plot_phi_init = DA.createVariable('truth_plot_phi_init','f8',('grid_pointS','one'))
  background_init = DA.createVariable('background_init','f8',('total data','number of members'))
  background_init_phi_c = DA.createVariable('background_init_phi_c','f8',('grid_pointS','twoS','number of members'))
  background_init_phi = DA.createVariable('background_init_phi','f8',('grid_pointS','number of members'))

  start_time = datetime.now()
  bg_error = {}                                               # background error dictionary
  an_error = {}                                               # analysis error dictionary
  for var in ['u','h','r']:
    bg_error[var] = np.zeros((ncyc))                   
    an_error[var] = np.zeros((ncyc))

  truth_plot_init[:,:], truth_plot_phi_c_init[:,:,:], truth_plot_phi_init[:,:] = Initializenew__(n,1,mindex,nindex)
  truth_plot[:,:,0], truth_plot_phi_c[:,:,:,0], truth_plot_phi[:,:,0]  = sweqssamp__(truth_plot_init[:,:] ,truth_plot_phi_init[:,:],truth_plot_phi_c_init[:,:,:],nsteps=dte-1,time_amp=0,n=n,n_ens=1,mindex=mindex,nindex=nindex,time=0,batch_number=0,members_per_batch=0)                # generate the truth model simulation. Time parameter = 0 to install mountain since first iteration
  truth = truth_plot[:,:,0]                # do we need these 3 lines? 
  truth_phi_c = truth_plot_phi_c[:,:,:,0]
  truth_phi = truth_plot_phi[:,:,0]
  
  background_init[:,:], background_init_phi_c[:,:,:], background_init_phi[:,:] = Initializenew__(n,nens,mindex,nindex)  
  background[:,:,0], background_phi_c[:,:,:,0], background_phi[:,:,0] = sweqssamp__(background_init[:,:],background_init_phi[:,:],background_init_phi_c[:,:,:],nsteps=dte-1,time_amp=0,n=n,n_ens=nens,mindex=mindex,nindex=nindex,time=0,batch_number=0,members_per_batch=0)             # propagate this ensemble. Time parameter = 0 to install mountain since forst iteration
  model = background[:,:,0]                # do we need these 3 lines? 
  model_phi_c = background_phi_c[:,:,:,0]
  model_phi = background_phi[:,:,0]

  for i in range(ncyc):                                       # go through ncyc DA cycles
      print('now doing cycle: '+str(i))
      for j,var in enumerate(['u','h','r']):                  # append background error 
         bg_error[var][i] = rmse(truth[j*n:n*(j+1),0],np.mean(model[j*n:n*(j+1),:],axis=1))
      truth_obs = truth.reshape(-1)                           # reshape makes shape of array (3*n) instead of (3*n,1). Needed for get_obs function
      obs_ncyc = i 
      obs = get_obs(truth=truth_obs,seed=obs_ncyc)            # create observations
      observation[i,:] = obs
      model = DA_method(obs,model,obs_ncyc)                   # compute analysis
      analysis[i,:,:] = model
     
      for j,var in enumerate(['u','h','r']):                  # append analysis error
          an_error[var][i] = rmse(truth[j*n:n*(j+1),0],np.mean(model[j*n:n*(j+1),:],axis=1))

      truth_plot[:,:,i+1], truth_plot_phi_c[:,:,:,i+1], truth_plot_phi[:,:,i+1] = sweqssamp__(truth,truth_phi, truth_phi_c,nsteps=dte,time_amp=(1+i)*dte,n=n,n_ens=1,mindex=mindex,nindex=nindex,time=(1+i)*dte,batch_number=0,members_per_batch=0)
      truth = truth_plot[:,:,i+1]
      truth_phi_c = truth_plot_phi_c[:,:,:,i+1]
      truth_phi = truth_plot_phi[:,:,i+1]

      background[:,:,i+1], background_phi_c[:,:,:,i+1], background_phi[:,:,i+1] = sweqssamp__(model, model_phi, model_phi_c,nsteps=dte,time_amp=(1+i)*dte,n=n,n_ens=nens,mindex=mindex,nindex=nindex,time=(1+i)*dte,batch_number=0,members_per_batch=0)
      model = background[:,:,i+1]
      model_phi_c = background_phi_c[:,:,:,i+1]
      model_phi = background_phi[:,:,i+1]

  pickle_obj = open(str(save_directory)+'dicts.an_error','wb')                  # saving an_error and bg_error in dictionary
  pickle.dump(an_error,pickle_obj)
  pickle_obj.close()
  pickle_obj = open(str(save_directory)+'dicts.bg_error','wb')
  pickle.dump(bg_error,pickle_obj)
  pickle_obj.close()
      
  time = datetime.now()-start_time
  print('time for DA initialisation with nens= ',nens,' and ncyc= ',ncyc,' equals= ', time)

  SWM.close()

  return
