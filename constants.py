#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 11:45:23 2019

@author: Yvonne.Ruckstuhl, edited by Kirsten 
"""
import numpy as np

n = 1000                                                    # number of gridpoints
dte = 75                                                    # time step between assimilations of observations in DA, corresponds to 5 model minutes
nens = 500                                                  # original number of ensemble members (number used for DA)                  
ncyc = 50                                                   # number of DA cycles
split = 100                                                 # used for separating ensemble so equal number(split) of members in each data file
multiply = 200                                              # factor by which to increase ensemble size after DA initialisation. Each data file will have split*multiply number of members
save_directory = ''                                         # where to save netCDF data files                                   
dx = 500.0                                                  # horizontal resolution of 500
dy = 500.0                                                  # vertical resolution
dts = 4.0                                                   # timestep used in discretisation
g = 10                                                      # gravitational constant
h_cloud = 90.02                                             # height threshold for cloud formation
phic = 899.77                                               # geopotential constant value above first threshold to allow for unstable convection.
h_0 = 90.0                                                  # resting height of fluid
gamma = g*h_0                                               # weight for negative buoyancy of rain (c in paper) 
h_rain = 90.4                                               # height threshold for rain
beta = 0.1                                                  # lag between cloud and rain formation.
alpha = 0.00014                                             # half-life of influence of rain of roughly 1 hour 
kr = 10.0                                                   # diffusion constant for rain.
hw = 5                                                      # half width of mountain ridge
amp = 1.2                                                   # amplitude of mountain ridge
mu = n/2                                                    # centre point of mountain ridge
t = 60                                                      # Number of timesteps in each netCDF file.
filter_parameter = 0.7                                      # RAW filter reduces computational mode.
alpha_filt = 0.53                                           # Want just above 0.5 to be conditionally stable. For RAW filter
tf = 60                                                     # timesteps between saved timesteps. This corresponds to 4 model minutes 


sig = [0.001, 0.01, 1.5]                                    # For DA. Noise added to truth to create observations
H = np.arange(3*n)
loc_radius = 4                                              # localisation radius for DA
mindex = 0                                                  # mountain index - 0: flat orography, 1: bell-shaped mountain, 2: power spectrum with sine envelope
nindex = 1                                                  # noise index - 0: no random pertubational noise in boundary layer, 1: random noise
n_array = np.array([[4,0.00895,0,n]])
#n_array = np.array([[4,0.005,0,n-1],[6,0.00,60,90]])    # info of random noise if nindex = 1. Each row indicates a location of noise in the order 
                                                         # [sigma (half width) of the noise field,amplitude, start, end of location where the model can randomnly selcet the noise position 
                                                         # Currently perturbations in whole domain like original (first row). 
                                                         # Standard amplitude = 0.00895. change to 0.000008 if mindex = 1.

rdiag = np.zeros((3*n),dtype='float32')                  # for DA
for i in range(3):
    rdiag[i*n:(i+1)*n] = sig[i]**2
rdiag = rdiag[H]
