#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 11:45:23 2019

@author: Yvonne.Ruckstuhl, edited by Kirsten ;) 
"""
import numpy as np
                                        
n = 1000                                                    # number of gridpoints
dte = 75                                                    # model time step for DA initialisation
nens = 500                                                  # original number of ensemble members                   
ncyc = 288                                                  # number of DA cycles = 24 hours
split = 100                                                 # number of members in each batch (before they are multiplied)
multiply = 10                                               # factor by which to increase ensemble size after DA initialisation
save_directory = '/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/'  # directory for saving output
dx = 500.0                                                  # horizontal resolution of 500
dy = 500.0                                                  # vertical resolution
dts = 4.0                                                   # timestep used in discretisation
g = 10                                                      # gravitational constant
h_cloud = 38.02                                             # height threshold for cloud formation
phic = 379.77                                               # geopotential constant value above first threshold to allow for unstable convection. Usually 899.77
h_0 = 38.0                                                  # resting height of fluid
gamma = g*h_0                                               # weight for negative buoyancy of rain (c in paper) originally 1.0
h_rain = 38.4                                               # height threshold for rain
beta = 0.1                                                  # lag between cloud and rain formation. standard 1/300
alpha = 0.00014                                             # half-life of influence of rain of roughly 1 hour was 0.00067
kr = 10.0                                                   # diffusion constant for rain. normally 50.0
hw = 5                                                      # half width of mountain ridge
amp = 1.2                                                   # amplitude of mountain ridge
mu = n/2                                                    # centre point of mountain ridge
t = 60 #360 # 720 # 1000    (720 is 48 hours)                   # length of the model simulation which uses discretisation timestep dts. Or timesteps in each netCDF file.
filter_parameter = 0.7                                      # usually 0.1. For RAW filter.
alpha_filt = 0.53                                           # usually 0.53. Want just above 0.5 to be conditionally stable. For RAW filter
tf = 60                                                     # timesteps between saved timesteps 
t_whole = 360 #360

# for new phi_c equation: 
#S_rad = 0 #-0.000015     #-0.0002
#S_for = 0
#S_for_2 = 0 #-0.00004 
gamma_2 = -2000         
freq = 1 # How many complete cycles of the diurnal cycle do I want in 24 hours? If 1 cycle, freq = 1.
for_loc = 0 #0.00001
time_for_gp = 14 # number of timesteps for a front (strong forcing) to cross1 dx. 
for_amp = 0 #0.00001
r_threshold = 0.000009

sig = [0.001, 0.01, 1.5]                                 # [0.001,  0.01, 0.0001]
H = np.arange(3*n)
loc_radius = 4                                              # normally 3
mindex = 0                                                  # mountain index - 0: flat orography, 1: bell-shaped mountain, 2: power spectrum with sine envelope
nindex = 1                                                  # noise index - 0: no random pertubational noise in boundary layer, 1: random noise
#n_array = np.array([[4,0.00895,0,n]])
n_array = np.array([[4,0,0,n]])    # amp normally 0.00895   To be used only during DA and init for edited run 4 and up. 
n_arrayamp = np.array([[4,0.011,0,n]]) # np.array([[4,0.00895,0,n]]) #  To be used only during free-run for edited run 4 and up
#n_array = np.array([[4,0.005,0,n-1],[6,0.00,60,90]])    # info of random noise if nindex = 1 Each row indicates a location of noise in the order 
                                                         # [sigma (half width) of the noise field,amplitude, start, end of place for model to choose 
                                                         # of noise field]. Currently perturbations in whole domain like original (first row). 
                                                         # Standard amplitude = 0.005. change to 0.000005 if mindex = 1.

rdiag = np.zeros((3*n),dtype='float32')
for i in range(3):
    rdiag[i*n:(i+1)*n] = sig[i]**2
rdiag = rdiag[H]
