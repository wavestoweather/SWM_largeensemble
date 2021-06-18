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
from DA_2019 import *
from scipy import stats
from datetime import datetime
from numpy import ndarray
from netCDF4 import Dataset
import numpy.ma as ma
import pickle
from datetime import datetime

# note those functions appended with a '__' save every single timestep. Those functions with 'ss__' save data every 'tf' timesteps.

def rmse(x,y):
  return (np.mean((x-y)**2))**0.5

def sweq__(x,nsteps,n,n_ens,mindex,nindex,time,batch_number,members_per_batch):    
    """
    x is the concatenated state of size 3*n. This is the relevant model input and ouput.
    nsteps is the number of timesteps to compute
    n is the total grid size
    n_ens is the number of ensemble members
    mindex = mountain index - 1: bell shaped mountain
    nindex = noise index - 1: stochastic Gaussian noise defined in constants
    time dicates time step count. Used to install orography and keep count for random seed
	  batch_number = corresponds to which file saving data in if using the splitting functionality to split the ensmeble up into different data files
	  members_per_batch = number of members to have in each of these 'batches' 
    """  
    phi = np.zeros((1,n+2,n_ens))                               # geopotential
    beta_new = np.zeros((n+2,n_ens))                            # array determines how much rain put in
    harray = np.zeros((n+2,n_ens))                              # orography
    y = np.zeros((3*n,n_ens,nsteps))                            # final outputting array
    ku = 2000                                                   # diffusion constant for u
    kh = 6000                                                   # diffusion coefficient for h
    
    if mindex == 1:                                             # bell shaped mountain ridge [Note: need to edit/check]
      ku = 3000                                                 # diffusion constant for u
      kh = ku                                                   # diffusion coefficient for depth h
      for i in range(n):
        harray[i+1,:] = amp*hw*hw/((float(i)-mu)*(float(i)-mu)+hw*hw)   # every ensemble member has same mountain

    if mindex == 2:                                             # k**(-1) power spectrum  [Note: need to edit/check]
      x_sin = np.arange(0,n)                                    # mindex #1 and #2 need work and checking
      s = 0.2+np.sin((x_sin-100)/500)
      np.random.seed(2)
      data = np.random.rand(n)-0.5                              # data centred around zero
      ps = np.abs(np.fft.fft(data,norm='ortho'))**2             # calculating the power spectrum
      for i in range(n_ens):
        harray[1:n+1,i]=ps*s

    harray[0,:] = harray[n,:]
    harray[n+1,:] = harray[1,:]

    u = np.zeros((3,n+2,n_ens))                                             # horizontal velocity 
    h = np.zeros((3,n+2,n_ens))                                             # depth 
    r = np.zeros((3,n+2,n_ens))                                             # rain
    u[0,1:n+1,:],u[1,1:n+1,:] = x[0:n,:],x[0:n,:]                           # filling in the middle of the u,h and r matrixes with input 
    h[0,1:n+1,:],h[1,1:n+1,:] = x[n:2*n,:],x[n:2*n,:]
    r[0,1:n+1,:],r[1,1:n+1,:] = x[2*n:3*n,:],x[2*n:3*n,:]
    
    u[1,0,:],u[0,0,:] = u[0,n,:],u[0,n,:]                                   # filling in the 2 outer columns of the u, h and r matrixes
    u[1,n+1,:],u[0,n+1,:] = u[0,1,:],u[0,1,:]                               # for boundary conditions from which the model effectively works. 
    h[1,0,:],h[0,0,:] = h[0,n,:],h[0,n,:]
    h[1,n+1,:],h[0,n+1,:] =  h[0,1,:],h[0,1,:]
    r[1,0,:],r[0,0,:] = r[0,n,:],r[0,n,:]
    r[1,n+1,:],r[0,n+1,:] = r[0,1,:],r[0,1,:]

    if time == 0:                                                           # intall orography - only happens at first time step
      h[:,:,:] = h[:,:,:] - harray[:,:]     

    for it in range(nsteps):                                                # loop over the given model timestep. Saves output after nsteps
      if nindex == 1:
        noise_seed = time+it
        u = noise(u,n,n_ens,n_array,noise_seed,batch_number,members_per_batch)                             # returns the u matrix with the 2nd row with stochastic noise across the ensemble
        u[1,0,:] = u[1,n,:]
        u[1,n+1,:] = u[1,1,:]

      if mindex == 1 or mindex==2:                                           # ignore mindex 
        phi[0,1:n+1,:] = np.where( h[1,1:n+1,:]+harray[1:n+1,:] > h_cloud , phic+g*harray[1:n+1,:], g*h[1,1:n+1,:] ) # if condition met, return phic in h matrix. If not, return g*h thing which would be the normal geopotential below the first threshold
      else:
        phi[0,1:n+1,:] = np.where( h[1,1:n+1,:] > h_cloud , phic, g*h[1,1:n+1,:] )

      phi[0,0,:]   = phi[0,n,:]                                             # boundary conditions for geopotential
      phi[0,n+1,:] = phi[0,1,:]
      phi[0,:,:] = phi[0,:,:] + gamma * r[1,:,:] 

      # shallow water equations
      u[2,1:n+1,:] = u[0,1:n+1,:] - (dts/(2*dx))*(u[1,2:n+2,:]**2 - u[1,0:n,:]**2) - (2*dts/dx)*(phi[0,1:n+1,:]-phi[0,0:n,:])  + (ku/(dx*dx))*(u[0,2:n+2,:] - 2*u[0,1:n+1,:] + u[0,0:n,:])*dts*2     # momentum equation 
      h[2,1:n+1,:] = h[0,1:n+1,:] - (dts/dx)*(u[1,2:n+2,:]*(h[1,1:n+1,:]+h[1,2:n+2,:]) - u[1,1:n+1,:]*(h[1,0:n,:]+h[1,1:n+1,:])) + (kh/(dx*dx))*(h[0,2:n+2,:] - 2*h[0,1:n+1,:] + h[0,0:n,:])*dts*2   # continuity equation 

      if mindex == 1 or mindex==2:
        mask = np.logical_and(h[1,1:n+1,:]+harray[1:n+1,:] > h_rain, u[1,2:n+2,:]-u[1,1:n+1,:] < 0)    # conditions for rain
      else:
        mask = np.logical_and(h[1,1:n+1,:] > h_rain, u[1,2:n+2,:]-u[1,1:n+1,:] < 0)   

      beta_new[1:n+1,:] = np.where( mask, beta , 0 )
      r[2,1:n+1,:] = r[0,1:n+1,:] - (dts/(2*dx))*(u[1,2:n+2,:]+u[1,1:n+1,:])*(r[1,2:n+2,:]-r[1,0:n,:]) - alpha*dts*2.0*r[1,1:n+1,:]-2.0*beta_new[1:n+1,:]*(dts/dx)*(u[1,2:n+2,:]-u[1,1:n+1,:]) + (kr/(dx*dx))*(r[0,2:n+2,:] - 2.0*r[0,1:n+1,:] + r[0,0:n,:])*dts*2 # rain equation  # with advection
      #r[2,1:n+1,:] = r[0,1:n+1,:] - alpha*dts*2.0*r[1,1:n+1,:]-2.0*beta_new[1:n+1,:]*(dts/dx)*(u[1,2:n+2,:]-u[1,1:n+1,:]) + (kr/(dx*dx))*(r[0,2:n+2,:] - 2.0*r[0,1:n+1,:] + r[0,0:n,:])*dts*2 # rain equation, no advection
      #r[2,1:n+1,:] = r[0,1:n+1,:] - alpha*dts*2.0*r[1,1:n+1,:]-2.0*beta_new[1:n+1,:]*(dts/dx)*(u[1,2:n+2,:]-u[1,1:n+1,:]) + (kr/(dx*dx))*(r[0,2:n+2,:] - 2.0*r[0,1:n+1,:] + r[0,0:n,:])*dts*2  - (dts/(dx))*(u[1,1:n+1,:])*(r[1,2:n+2,:]-r[1,0:n,:])  # other advection term 
     
      u[2,0,:]   = u[2,n,:]                                                  # boundary conditions
      u[2,n+1,:] = u[2,1,:]
      h[2,0,:]   = h[2,n,:]
      h[2,n+1,:] = h[2,1,:]
      r[2,0,:]   = r[2,n,:]
      r[2,n+1,:] = r[2,1,:]
      r[2,:,:] = np.where(r[2,:,:]<0.,0,r[2,:,:])

      d = filter_parameter*.5*(u[2,:,:] - 2.*u[1,:,:] + u[0,:,:])            # RAW filter. Accounts for the growing computational mode.
      u[0,:,:] = u[1,:,:] + alpha_filt*d
      u[1,:,:] = u[2,:,:] - (1-alpha_filt)*d                     
      d = filter_parameter*.5*(h[2,:,:] - 2.*h[1,:,:] + h[0,:,:])
      h[0,:,:] = h[1,:,:] + alpha_filt*d
      h[1,:,:] = h[2,:,:] - (1-alpha_filt)*d
      d = filter_parameter*.5*(r[2,:,:] - 2.*r[1,:,:] + r[0,:,:])
      r[0,:,:] = r[1,:,:] + alpha_filt*d
      r[1,:,:] = r[2,:,:] - (1-alpha_filt)*d

      y[0:n,:,it] = u[2,1:n+1,:]
      y[n:2*n,:,it] = h[2,1:n+1,:]
      y[2*n:3*n,:,it] = r[2,1:n+1,:] 

    return y

def sweqss__(x,nsteps,n,n_ens,mindex,nindex,time,batch_number,members_per_batch): # ss for selected saving (otherwise same as sweq__ function)
    """
    x is the concatenated state of size 3*n. This is the relevant model input and ouput.
    nsteps is the number of timesteps
    n is the total grid size
    n_ens is the number of ensemble members
    mindex = mountain index - 1: bell shaped mountain
    nindex = noise index - 1: stochastic Gaussian noise throughout grid domain
    time indicates whether on first time step. Used to install orography and keep count for random seed
	  batch_number = corresponds to which file saving data in if have used the splittng functionality to split the ensmeble up into different data files
	  members_per_batch = number of members to have in each of these 'batches'     
	"""  
    phi = np.zeros((1,n+2,n_ens))                               # geopotential
    beta_new = np.zeros((n+2,n_ens))                            # controls amount of rain added to model
    harray = np.zeros((n+2,n_ens))                              # orography
    y = np.zeros((3*n,n_ens))                                   # final outputting array
    ku = 2000                                                   # diffusion constant for u usually 3000
    kh = 6000                                                   # diffusion coefficient for depth h
    
    if mindex == 1:                                             # bell shaped mountain ridge.   # ignore for now (need to edit)
      ku = 3000                                                 # diffusion constant for u
      kh = ku                                                   # diffusion coefficient for depth h
      for i in range(n):
        harray[i+1,:] = amp*hw*hw/((float(i)-mu)*(float(i)-mu)+hw*hw)   # every ensemble member has same mountain

    if mindex == 2:                                             # k**(-1) power spectrum.        # ignore for now (need to edit)
      x_sin = np.arange(0,n)                                    # mindex #1 and #2 need work and checking
      s = 0.2+np.sin((x_sin-100)/500)
      np.random.seed(2)
      data = np.random.rand(n)-0.5                              # data centred around zero
      ps = np.abs(np.fft.fft(data,norm='ortho'))**2             # calculating the power spectrum
      for i in range(n_ens):
        harray[1:n+1,i]=ps*s

    harray[0,:] = harray[n,:]
    harray[n+1,:] = harray[1,:]

    u = np.zeros((3,n+2,n_ens))                                             # horizontal velocity 
    h = np.zeros((3,n+2,n_ens))                                             # depth 
    r = np.zeros((3,n+2,n_ens))                                             # rain 
    u[0,1:n+1,:],u[1,1:n+1,:] = x[0:n,:],x[0:n,:]                           # filling in the middle of the u,h and r matrixes with input 
    h[0,1:n+1,:],h[1,1:n+1,:] = x[n:2*n,:],x[n:2*n,:]
    r[0,1:n+1,:],r[1,1:n+1,:] = x[2*n:3*n,:],x[2*n:3*n,:]
    
    u[1,0,:],u[0,0,:] = u[0,n,:],u[0,n,:]                                   # filling in the 2 outer columns of the u, h and r matrixes
    u[1,n+1,:],u[0,n+1,:] = u[0,1,:],u[0,1,:]                               # for boundary conditions from which the model effectively works. 
    h[1,0,:],h[0,0,:] = h[0,n,:],h[0,n,:]
    h[1,n+1,:],h[0,n+1,:] =  h[0,1,:],h[0,1,:]
    r[1,0,:],r[0,0,:] = r[0,n,:],r[0,n,:]
    r[1,n+1,:],r[0,n+1,:] = r[0,1,:],r[0,1,:]

    if time == 0:                                                           # intall orography - only happens at first time step
      h[:,:,:] = h[:,:,:] - harray[:,:]     

    for it in range(nsteps):                                                # loop over the given model timestep. Saves output after nsteps
      if nindex == 1:
        noise_seed = time+it
        u = noise(u,n,n_ens,n_array,noise_seed,batch_number,members_per_batch)                             # returns the u matrix with the 2nd row with stochastic noise across the ensemble
        u[1,0,:] = u[1,n,:]
        u[1,n+1,:] = u[1,1,:]

      if mindex == 1 or mindex==2:
        phi[0,1:n+1,:] = np.where( h[1,1:n+1,:]+harray[1:n+1,:] > h_cloud , phic+g*harray[1:n+1,:], g*h[1,1:n+1,:] ) # if condition met, return phic in h matrix. If not, return g*h thing which would be the normal geopotential below the first threshold
      else:
        phi[0,1:n+1,:] = np.where( h[1,1:n+1,:] > h_cloud , phic, g*h[1,1:n+1,:] )

      phi[0,0,:]   = phi[0,n,:]                                             # boundary conditions for geopotential
      phi[0,n+1,:] = phi[0,1,:]
      phi[0,:,:] = phi[0,:,:] + gamma * r[1,:,:] 

      # shallow water equations 
      u[2,1:n+1,:] = u[0,1:n+1,:] - (dts/(2*dx))*(u[1,2:n+2,:]**2 - u[1,0:n,:]**2) - (2*dts/dx)*(phi[0,1:n+1,:]-phi[0,0:n,:])  + (ku/(dx*dx))*(u[0,2:n+2,:] - 2*u[0,1:n+1,:] + u[0,0:n,:])*dts*2     # momentum equation
      h[2,1:n+1,:] = h[0,1:n+1,:] - (dts/dx)*(u[1,2:n+2,:]*(h[1,1:n+1,:]+h[1,2:n+2,:]) - u[1,1:n+1,:]*(h[1,0:n,:]+h[1,1:n+1,:])) + (kh/(dx*dx))*(h[0,2:n+2,:] - 2*h[0,1:n+1,:] + h[0,0:n,:])*dts*2   # continuity equation 

      if mindex == 1 or mindex==2:
        mask = np.logical_and(h[1,1:n+1,:]+harray[1:n+1,:] > h_rain, u[1,2:n+2,:]-u[1,1:n+1,:] < 0)    # conditions for rain
      else:
        mask = np.logical_and(h[1,1:n+1,:] > h_rain, u[1,2:n+2,:]-u[1,1:n+1,:] < 0)   

      beta_new[1:n+1,:] = np.where( mask, beta , 0 )
      r[2,1:n+1,:] = r[0,1:n+1,:] - (dts/(2*dx))*(u[1,2:n+2,:]+u[1,1:n+1,:])*(r[1,2:n+2,:]-r[1,0:n,:]) - alpha*dts*2.0*r[1,1:n+1,:]-2.0*beta_new[1:n+1,:]*(dts/dx)*(u[1,2:n+2,:]-u[1,1:n+1,:]) + (kr/(dx*dx))*(r[0,2:n+2,:] - 2.0*r[0,1:n+1,:] + r[0,0:n,:])*dts*2 # rain equation  # with advection
      #r[2,1:n+1,:] = r[0,1:n+1,:] - alpha*dts*2.0*r[1,1:n+1,:]-2.0*beta_new[1:n+1,:]*(dts/dx)*(u[1,2:n+2,:]-u[1,1:n+1,:]) + (kr/(dx*dx))*(r[0,2:n+2,:] - 2.0*r[0,1:n+1,:] + r[0,0:n,:])*dts*2 # rain equation, no advection
      #r[2,1:n+1,:] = r[0,1:n+1,:] - alpha*dts*2.0*r[1,1:n+1,:]-2.0*beta_new[1:n+1,:]*(dts/dx)*(u[1,2:n+2,:]-u[1,1:n+1,:]) + (kr/(dx*dx))*(r[0,2:n+2,:] - 2.0*r[0,1:n+1,:] + r[0,0:n,:])*dts*2  - (dts/(dx))*(u[1,1:n+1,:])*(r[1,2:n+2,:]-r[1,0:n,:])  # other advection term 
     
      u[2,0,:]   = u[2,n,:]                                                  # boundary conditions
      u[2,n+1,:] = u[2,1,:]
      h[2,0,:]   = h[2,n,:]
      h[2,n+1,:] = h[2,1,:]
      r[2,0,:]   = r[2,n,:]
      r[2,n+1,:] = r[2,1,:]
      r[2,:,:] = np.where(r[2,:,:]<0.,0,r[2,:,:])

      d = filter_parameter*.5*(u[2,:,:] - 2.*u[1,:,:] + u[0,:,:])            # RAW filter. Accounts for the growing computational mode.
      u[0,:,:] = u[1,:,:] + alpha_filt*d
      u[1,:,:] = u[2,:,:] - (1-alpha_filt)*d                     
      d = filter_parameter*.5*(h[2,:,:] - 2.*h[1,:,:] + h[0,:,:])
      h[0,:,:] = h[1,:,:] + alpha_filt*d
      h[1,:,:] = h[2,:,:] - (1-alpha_filt)*d
      d = filter_parameter*.5*(r[2,:,:] - 2.*r[1,:,:] + r[0,:,:])
      r[0,:,:] = r[1,:,:] + alpha_filt*d
      r[1,:,:] = r[2,:,:] - (1-alpha_filt)*d

    y[0:n,:] = u[2,1:n+1,:]
    y[n:2*n,:] = h[2,1:n+1,:]
    y[2*n:3*n,:] = r[2,1:n+1,:] 

    return y


def Initialize__(n,nens,mindex,nindex):
  '''
  Initialise x array input for sweq functions when initialising DA. 
  n - grid size 
  nens - number of ensembles. Note!: no multiplying factor here since 'Initialize' used for DA only.
  mindex - orography added to model (set to zero)
  nindex - random perturbations added to model
  '''
  state = np.zeros((3*n,nens))
  state[0:n] = np.zeros((n,nens))+0.                          # velocity 0m/s
  state[n:2*n] = np.zeros((n,nens))+90.                       # fluid depth 90m

  state = sweq__(state,1000,n,nens,mindex,nindex,time=0,batch_number=0,members_per_batch=0)      # state to begin with. 1000 nsteps taken
  return state

def noise(u,n,n_ens,n_array,noise_seed,batch_number,members_per_batch):   
    ''' 
    Random stochastic noise in form of convergence added to the velocity field. Mimics turbulence in the boundary layer. 
    n_array - tells model where noise is located. narray.shape = (number of perturbation areas, features of pert area)
    features of pert area = [half width, amplitude, start, end of noise field region]
    seeded so 'random' noise different on each iteration cycle and every ensemble member
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
        e_update = e + (batch_number*members_per_batch)               # each batch doesn't recieve the same forcing
        noise_seed_chosen = ((noise_seed+e_update)*(e_update+2))-(noise_seed+e_update*3)   # choose the random seed so that each ensemble member is perturbed differently, and at different iterations differently also. 
        while noise_seed_chosen>=2**32:
            noise_seed_chosen = noise_seed_chosen-(2**32)             # returns error if seed above this value
        np.random.seed(noise_seed_chosen)  
        pos = np.random.randint(start,end)                            # each ensemble will have own random noise added within bounds created
        unoise[pos:pos+n,e] = unoise[pos:pos+n,e] + zsum
        unoise[0:n,e] = unoise[0:n,e] + unoise[n:n+n,e]
        unoise[n:n+n,e] = 0 
        u[1,1:n+1,e] = u[1,1:n+1,e] + unoise[0:n,e]
    return u

def model_init__(start_time,batch_number,members_per_batch):
  '''
  Running the SWE with the IC from DA_init_vis_dt__
  start_time - tracks time (for seed)
  batch_number - as in sweq__
  members_per_batch - as in sweq__
  Creates 'actual_truth' run which is the truth from the DA initialisation continued. So it really has a total ensemble size of (nens*multiply)+1
  'Truth' is the ensemble member runs
  If multiply = True, then enemble multiplied. Don't want this if continuing run! 
  Need to change actual_truth and DA_time dependent on which DA used. 
  '''
  start_time_model = datetime.now()
  if start_time==0:
      SWMDA = Dataset(str(save_directory)+'SWMDA.nc','r')                # read DA netCDF.  # Change this depending on own labelling of files
      SWM = Dataset(str(save_directory)+'SWM0_'+str(batch_number)+'_'+str(members_per_batch)+'.nc','w')    # create a new file, don't use same file as DA anymore!!!
  else:
      OLD = Dataset(str(save_directory)+'SWM'+str(start_time-t)+'_'+str(batch_number)+'_'+str(members_per_batch)+'.nc','r')
      SWM = Dataset(str(save_directory)+'SWM'+str(start_time)+'_'+str(batch_number)+'_'+str(members_per_batch)+'.nc','w',format='NETCDF4') 
      
  free_run = SWM.createGroup('free_run')                                   # create free_run group
  no_ens = free_run.createDimension('number of members',members_per_batch) # create dimensions
  time = free_run.createDimension('time',None)
  total_data = free_run.createDimension('total data',3*n)
  one = free_run.createDimension('one',1)
  truth = free_run.createVariable('truth','f8',('total data','number of members','time')) # create variables # 'truth' - ensemble model members
  actual_truth = free_run.createVariable('actual_truth','f8',('total data','one','time')) # 'actual_truth' = truth

  DA_time = (1+ncyc)*(dte)+1000                                                          # keep track of time for seed
  if start_time==0:             
      actual_truth[:,:,0] = ma.getdata(SWMDA['/DA/truth_plot'][:,:,(ncyc*dte)-1])        # picking up the actual_truth from the DA initialisation
  else:
      actual_truth[:,:,0] = ma.getdata(OLD['/free_run/actual_truth'][:,:,t-1])

  actual_truth[:,:,1:] = sweq__(actual_truth[:,:,0],nsteps=t-1,n=n,n_ens=1,mindex=mindex,nindex=nindex,time=(start_time)+DA_time,batch_number=0,members_per_batch=0) # integrating this 'real' truth further in time for model simulation

  if start_time ==0:                                                                                                          # first time step of truth given by DA initialisation
      truth[:,:,0] = np.repeat(ma.getdata(SWMDA['/DA/analysis'][ncyc-1,:,batch_number*split:(batch_number+1)*split]),repeats=multiply,axis=1)          # change this line depending on how splitting up analysis!                 # method depends on whether or not you wish to expand the ensemble size
      if len(truth[1,:,1]) == members_per_batch:
      	print('split correctly')
      else:
      	print('not split correctly!!!!')
  else:
      truth[:,:,0] = ma.getdata(OLD['/free_run/truth'][:,:,t-1])

  truth[:,:,1:] = sweq__(truth[:,:,0],nsteps=t-1,n=n,n_ens=members_per_batch,mindex=mindex,nindex=nindex,time=start_time+DA_time,batch_number=batch_number,members_per_batch=members_per_batch)     # generate model simulation, store in 'truth'

  time = datetime.now()-start_time_model
  print('time for model run with nens= ',members_per_batch,' equals= ', time)
  return 

def model_initss__(start_time,batch_number,members_per_batch):   # ss for selected saving (otherwise same as model_init__ function)
  '''
  Running the SWE with the IC from DA_init_vis_dt__
  start_time - tracks time (for seed)
  batch_number - as in sweq__
  members_per_batch - as in sweq__
  Creates 'actual_truth' run which is the truth from the DA initialisation continued. So it really has a total ensemble size of (nens*multiply)+1
  'Truth' is the ensemble member runs
  If multiply = True, then enemble multiplied. Don't want this if continuing run! 
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
  truth = free_run.createVariable('truth','f8',('total data','number of members','time')) # create variables
  actual_truth = free_run.createVariable('actual_truth','f8',('total data','one','time'))

  DA_time = (1+ncyc)*(dte)+1000         
  if start_time==0:             
      actual_truth[:,:,0] = ma.getdata(SWMDA['/DA/truth_plot'][:,:,(ncyc*dte)-1])        # picking up the truth from the DA initialisation
  else:
      actual_truth[:,:,0] = ma.getdata(OLD['/free_run/actual_truth'][:,:,t-1])

  for i in range(t):
    print('Computing the '+str(i)+' of t now.')
    if i == 0:
      actual_truth[:,:,0] = sweqss__(actual_truth[:,:,0],nsteps=tf,n=n,n_ens=1,mindex=mindex,nindex=nindex,time=(start_time*tf)+DA_time,batch_number=0,members_per_batch=0) 
    else:
      actual_truth[:,:,i] = sweqss__(actual_truth[:,:,i-1],nsteps=tf,n=n,n_ens=1,mindex=mindex,nindex=nindex,time=(start_time*tf)+DA_time+(i*tf),batch_number=0,members_per_batch=0) # integrating this 'real' truth further in time for model simulation

  if start_time ==0:                                                                                                          # first time step of truth given by DA initialisation
      truth[:,:,0] = np.repeat(ma.getdata(SWMDA['/DA/analysis'][ncyc-1,:,batch_number*split:(batch_number+1)*split]),repeats=multiply,axis=1)          # change this line depending on how splitting up analysis!                 # method depends on whether or not you wish to expand the ensemble size
      if len(truth[1,:,1]) == members_per_batch:
        print('split correctly')
      else:
        print('not split correctly!!!!')
  else:
      truth[:,:,0] = ma.getdata(OLD['/free_run/truth'][:,:,t-1])

  for i in range(t):
    if i == 0:
      truth[:,:,0] = sweqss__(truth[:,:,0],nsteps=tf,n=n,n_ens=members_per_batch,mindex=mindex,nindex=nindex,time=(start_time*tf)+DA_time,batch_number=batch_number,members_per_batch=members_per_batch) 
    else:
      truth[:,:,i] = sweqss__(truth[:,:,i-1],nsteps=tf,n=n,n_ens=members_per_batch,mindex=mindex,nindex=nindex,time=(start_time*tf)+DA_time+(i*tf),batch_number=batch_number,members_per_batch=members_per_batch)     # generate model simulation, store in 'truth'

  time = datetime.now()-start_time_model
  print('time for model run with nens= ',members_per_batch,' equals= ', time)
  return 

def DA_init_vis_dt__(DA_method):
  '''
  Initialising the model with DA using the DA_method(usually EnKF). Goes through ncyc cycles before outputting the array model_start which are the IC. 
  Also returns rmse error of the background and analysis, bk_error and an_error
  Tracks time for initialisation so seeding
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
  truth_data_points = DA.createDimension('truth data points',dte*(ncyc+1))
  init_length = DA.createDimension('init length',1000)

  analysis = DA.createVariable('analysis','f8',('number of DA cycles','total data','number of members'))
  background = DA.createVariable('background','f8',('total data','number of members','truth data points'))
  observation = DA.createVariable('observation','f8',('number of DA cycles','total data'))
  truth_plot = DA.createVariable('truth_plot','f8',('total data','one','truth data points'))
  truth_plot_init = DA.createVariable('truth_plot_init','f8',('total data','one','init length'))
  background_init = DA.createVariable('background_init','f8',('total data','number of members','init length'))

  start_time = datetime.now()
  bg_error = {}                                               # background error dictionary
  an_error = {}                                               # analysis error dictionary
  for var in ['u','h','r']:
    bg_error[var] = np.zeros((ncyc))                   
    an_error[var] = np.zeros((ncyc))
   
  truth_plot_init[:,:,:] = Initialize__(n,1,mindex,nindex)
  truth_plot[:,:,0] = truth_plot_init[:,:,999]                # use last time step # truth run
  truth_plot[:,:,1:dte] = sweq__(truth_plot[:,:,0],nsteps=dte-1,n=n,n_ens=1,mindex=mindex,nindex=nindex,time=1000,batch_number=0,members_per_batch=0)                # generate the truth run. Time parameter = 0 to install mountain since first iteration (don't use for now)
  truth = truth_plot[:,:,dte-1]
  
  background_init[:,:,:] = Initialize__(n,nens,mindex,nindex) # background
  background[:,:,0] = background_init[:,:,999]                # generate the initial ensemble, use last time step
  background[:,:,1:dte] = sweq__(background[:,:,0],nsteps=dte-1,n=n,n_ens=nens,mindex=mindex,nindex=nindex,time=1000,batch_number=0,members_per_batch=0)             # propagate this ensemble. Time parameter = 0 to install mountain since forst iteration
  model = background[:,:,dte-1]

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

      truth_plot[:,:,((i*dte)+dte):((i*dte)+(2*dte))] = sweq__(truth,nsteps=dte,n=n,n_ens=1,mindex=mindex,nindex=nindex,time=1000+((1+i)*dte),batch_number=0,members_per_batch=0)
      truth = truth_plot[:,:,(2*dte)+(i*dte)-1]

      background[:,:,((i*dte)+dte):((i*dte)+(2*dte))] = sweq__(model,nsteps=dte,n=n,n_ens=nens,mindex=mindex,nindex=nindex,time=1000+((1+i)*dte),batch_number=0,members_per_batch=0)
      model = background[:,:,(2*dte)+(i*dte)-1]

  pickle_obj = open(str(save_directory)+'dicts.an_error','wb')                  # saving an_error and bg_error in dictionary
  pickle.dump(an_error,pickle_obj)
  pickle_obj.close()
  pickle_obj = open(str(save_directory)+'dicts.bg_error','wb')
  pickle.dump(bg_error,pickle_obj)
  pickle_obj.close()
      
  time = datetime.now()-start_time
  print('time for DA initialisation with nens= ',nens,' and ncyc= ',ncyc,' equals= ', time)

  return

