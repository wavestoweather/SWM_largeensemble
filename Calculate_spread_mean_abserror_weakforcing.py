import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse
import numpy.ma as ma
from netCDF4 import Dataset
import pickle

################################################

# calculate stats from the distributions. Do for individual grid points and then average toegther. 
# 1) mean
# 2) spread (stdev wrt mean of ensemble memebers)
# 3) absolute error of mean to truth

################################################

# import data
SWM00 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM0_0_1000.nc','r')
SWM01 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM0_1_1000.nc','r')
SWM02 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM0_2_1000.nc','r')
SWM03 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM0_3_1000.nc','r')
SWM04 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM0_4_1000.nc','r')

SWM600 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM60_0_1000.nc','r')
SWM601 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM60_1_1000.nc','r')
SWM602 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM60_2_1000.nc','r')
SWM603 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM60_3_1000.nc','r')
SWM604 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM60_4_1000.nc','r')

SWM1200 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM120_0_1000.nc','r')
SWM1201 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM120_1_1000.nc','r')
SWM1202 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM120_2_1000.nc','r')
SWM1203 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM120_3_1000.nc','r')
SWM1204 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM120_4_1000.nc','r')

SWM1800 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM180_0_1000.nc','r')
SWM1801 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM180_1_1000.nc','r')
SWM1802 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM180_2_1000.nc','r')
SWM1803 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM180_3_1000.nc','r')
SWM1804 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM180_4_1000.nc','r')

SWM2400 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM240_0_1000.nc','r')
SWM2401 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM240_1_1000.nc','r')
SWM2402 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM240_2_1000.nc','r')
SWM2403 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM240_3_1000.nc','r')
SWM2404 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM240_4_1000.nc','r')

SWM3000 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM300_0_1000.nc','r')
SWM3001 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM300_1_1000.nc','r')
SWM3002 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM300_2_1000.nc','r')
SWM3003 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM300_3_1000.nc','r')
SWM3004 = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM300_4_1000.nc','r')

dp = np.arange(0,1000)

    #########################################################

data_phi_c_0 = np.concatenate((ma.getdata(SWM00['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM600['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM1200['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM1800['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM2400['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM3000['/free_run/truth_phi_c'][dp,:,:,:])),axis=3)
data_phi_c_1 = np.concatenate((ma.getdata(SWM01['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM601['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM1201['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM1801['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM2401['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM3001['/free_run/truth_phi_c'][dp,:,:,:])),axis=3)
data_phi_c_2 = np.concatenate((ma.getdata(SWM02['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM602['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM1202['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM1802['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM2402['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM3002['/free_run/truth_phi_c'][dp,:,:,:])),axis=3)
data_phi_c_3 = np.concatenate((ma.getdata(SWM03['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM603['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM1203['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM1803['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM2403['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM3003['/free_run/truth_phi_c'][dp,:,:,:])),axis=3)
data_phi_c_4 = np.concatenate((ma.getdata(SWM04['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM604['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM1204['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM1804['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM2404['/free_run/truth_phi_c'][dp,:,:,:]),ma.getdata(SWM3004['/free_run/truth_phi_c'][dp,:,:,:])),axis=3)

data_phi_c = np.concatenate((data_phi_c_0,data_phi_c_1,data_phi_c_2,data_phi_c_3,data_phi_c_4),axis=2)

data_phi_c = data_phi_c[:,1,:,:]
# get the data to have dimension n
# phi data 

data_phi_0 = np.concatenate((ma.getdata(SWM00['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM600['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM1200['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM1800['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM2400['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM3000['/free_run/truth_phi'][dp,:,:])),axis=2)
data_phi_1 = np.concatenate((ma.getdata(SWM01['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM601['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM1201['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM1801['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM2401['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM3001['/free_run/truth_phi'][dp,:,:])),axis=2)
data_phi_2 = np.concatenate((ma.getdata(SWM02['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM602['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM1202['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM1802['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM2402['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM3002['/free_run/truth_phi'][dp,:,:])),axis=2)
data_phi_3 = np.concatenate((ma.getdata(SWM03['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM603['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM1203['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM1803['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM2403['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM3003['/free_run/truth_phi'][dp,:,:])),axis=2)
data_phi_4 = np.concatenate((ma.getdata(SWM04['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM604['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM1204['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM1804['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM2404['/free_run/truth_phi'][dp,:,:]),ma.getdata(SWM3004['/free_run/truth_phi'][dp,:,:])),axis=2)

data_phi = np.concatenate((data_phi_0,data_phi_1,data_phi_2,data_phi_3,data_phi_4),axis=1)

# wind data 

data_w_0 = np.concatenate((ma.getdata(SWM00['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM600['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1200['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1800['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM2400['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM3000['/free_run/truth'][0:1000,:,:])),axis=2)
data_w_1 = np.concatenate((ma.getdata(SWM01['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM601['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1201['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1801['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM2401['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM3001['/free_run/truth'][0:1000,:,:])),axis=2)
data_w_2 = np.concatenate((ma.getdata(SWM02['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM602['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1202['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1802['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM2402['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM3002['/free_run/truth'][0:1000,:,:])),axis=2)
data_w_3 = np.concatenate((ma.getdata(SWM03['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM603['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1203['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1803['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM2403['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM3003['/free_run/truth'][0:1000,:,:])),axis=2)
data_w_4 = np.concatenate((ma.getdata(SWM04['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM604['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1204['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1804['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM2404['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM3004['/free_run/truth'][0:1000,:,:])),axis=2)

data_w = np.concatenate((data_w_0,data_w_1,data_w_2,data_w_3,data_w_4),axis=1)

# interpolate wind data
n = 1000
a = np.insert(data_w,n,data_w[0,:,:],axis=0)
data_wind = 0.5*(np.sum([a[1:n+1,:,:]+a[0:n,:,:]],axis=0))

# height data

data_h_0 = np.concatenate((ma.getdata(SWM00['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM600['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM1200['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM1800['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM2400['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM3000['/free_run/truth'][dp + 1000,:,:])),axis=2)
data_h_1 = np.concatenate((ma.getdata(SWM01['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM601['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM1201['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM1801['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM2401['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM3001['/free_run/truth'][dp + 1000,:,:])),axis=2)
data_h_2 = np.concatenate((ma.getdata(SWM02['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM602['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM1202['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM1802['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM2402['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM3002['/free_run/truth'][dp + 1000,:,:])),axis=2)
data_h_3 = np.concatenate((ma.getdata(SWM03['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM603['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM1203['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM1803['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM2403['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM3003['/free_run/truth'][dp + 1000,:,:])),axis=2)
data_h_4 = np.concatenate((ma.getdata(SWM04['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM604['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM1204['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM1804['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM2404['/free_run/truth'][dp + 1000,:,:]),ma.getdata(SWM3004['/free_run/truth'][dp + 1000,:,:])),axis=2)

data_height = np.concatenate((data_h_0,data_h_1,data_h_2,data_h_3,data_h_4),axis=1)

# rain data

data_r_0 = np.concatenate((ma.getdata(SWM00['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM600['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM1200['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM1800['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM2400['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM3000['/free_run/truth'][dp + 2000,:,:])),axis=2)
data_r_1 = np.concatenate((ma.getdata(SWM01['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM601['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM1201['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM1801['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM2401['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM3001['/free_run/truth'][dp + 2000,:,:])),axis=2)
data_r_2 = np.concatenate((ma.getdata(SWM02['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM602['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM1202['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM1802['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM2402['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM3002['/free_run/truth'][dp + 2000,:,:])),axis=2)
data_r_3 = np.concatenate((ma.getdata(SWM03['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM603['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM1203['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM1803['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM2403['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM3003['/free_run/truth'][dp + 2000,:,:])),axis=2)
data_r_4 = np.concatenate((ma.getdata(SWM04['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM604['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM1204['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM1804['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM2404['/free_run/truth'][dp + 2000,:,:]),ma.getdata(SWM3004['/free_run/truth'][dp + 2000,:,:])),axis=2)

data_rain = np.concatenate((data_r_0,data_r_1,data_r_2,data_r_3,data_r_4),axis=1)

# truth

data_phi_c_truth = np.concatenate((ma.getdata(SWM00['/free_run/actual_truth_phi_c'][dp,:,:,:]),ma.getdata(SWM600['/free_run/actual_truth_phi_c'][dp,:,:,:]),ma.getdata(SWM1200['/free_run/actual_truth_phi_c'][dp,:,:,:]),ma.getdata(SWM1800['/free_run/actual_truth_phi_c'][dp,:,:,:]),ma.getdata(SWM2400['/free_run/actual_truth_phi_c'][dp,:,:,:]),ma.getdata(SWM3000['/free_run/actual_truth_phi_c'][dp,:,:,:])),axis=3)
data_phi_c_truth = data_phi_c_truth[:,1,:,:]
data_phi_truth = np.concatenate((ma.getdata(SWM00['/free_run/actual_truth_phi'][dp,:,:]),ma.getdata(SWM600['/free_run/actual_truth_phi'][dp,:,:]),ma.getdata(SWM1200['/free_run/actual_truth_phi'][dp,:,:]),ma.getdata(SWM1800['/free_run/actual_truth_phi'][dp,:,:]),ma.getdata(SWM2400['/free_run/actual_truth_phi'][dp,:,:]),ma.getdata(SWM3000['/free_run/actual_truth_phi'][dp,:,:])),axis=2)
data_wind_truth = np.concatenate((ma.getdata(SWM00['/free_run/actual_truth'][0:1000,:,:]),ma.getdata(SWM600['/free_run/actual_truth'][0:1000,:,:]),ma.getdata(SWM1200['/free_run/actual_truth'][0:1000,:,:]),ma.getdata(SWM1800['/free_run/actual_truth'][0:1000,:,:]),ma.getdata(SWM2400['/free_run/actual_truth'][0:1000,:,:]),ma.getdata(SWM3000['/free_run/actual_truth'][0:1000,:,:])),axis=2)
# interpolate wind data
n = 1000
a = np.insert(data_wind_truth,n,data_wind_truth[0,:,:],axis=0)
data_wind_truth = 0.5*(np.sum([a[1:n+1,:,:]+a[0:n,:,:]],axis=0))
data_height_truth = np.concatenate((ma.getdata(SWM00['/free_run/actual_truth'][dp + 1000,:,:]),ma.getdata(SWM600['/free_run/actual_truth'][dp + 1000,:,:]),ma.getdata(SWM1200['/free_run/actual_truth'][dp + 1000,:,:]),ma.getdata(SWM1800['/free_run/actual_truth'][dp + 1000,:,:]),ma.getdata(SWM2400['/free_run/actual_truth'][dp + 1000,:,:]),ma.getdata(SWM3000['/free_run/actual_truth'][dp + 1000,:,:])),axis=2)
data_rain_truth = np.concatenate((ma.getdata(SWM00['/free_run/actual_truth'][dp + 2000,:,:]),ma.getdata(SWM600['/free_run/actual_truth'][dp + 2000,:,:]),ma.getdata(SWM1200['/free_run/actual_truth'][dp + 2000,:,:]),ma.getdata(SWM1800['/free_run/actual_truth'][dp + 2000,:,:]),ma.getdata(SWM2400['/free_run/actual_truth'][dp + 2000,:,:]),ma.getdata(SWM3000['/free_run/actual_truth'][dp + 2000,:,:])),axis=2)

#########################################################
# create .nc file to put data in:
DIST = Dataset('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/stats/stats_weakforcing27.nc','w',format='NETCDF4')
# create groups
G_WIND = DIST.createGroup('G_WIND')
G_HEIGHT = DIST.createGroup('G_HEIGHT')
G_RAIN = DIST.createGroup('G_RAIN')
G_EXTRA = DIST.createGroup('G_EXTRA')
# create dimensions
#wind
no_ens = G_WIND.createDimension('stats',6)
time = G_WIND.createDimension('number of timesteps',360)
# height
no_ens = G_HEIGHT.createDimension('stats',6)
time = G_HEIGHT.createDimension('number of timesteps',360)
# rain
no_ens = G_RAIN.createDimension('stats',6)
time = G_RAIN.createDimension('number of timesteps',360)
# extra
no_ens = G_EXTRA.createDimension('stats',6)
time = G_EXTRA.createDimension('number of timesteps',360)

# create variables
# phi_c 
phi_c_stats = G_EXTRA.createVariable('phi_c_stats','f8',('stats','number of timesteps'))
# phi
phi_stats = G_EXTRA.createVariable('phi_stats','f8',('stats','number of timesteps'))
# wind
wind_stats = G_WIND.createVariable('wind_stats','f8',('stats','number of timesteps'))
# height
height_stats = G_HEIGHT.createVariable('height_stats','f8',('stats','number of timesteps'))
# rain
rain_stats = G_RAIN.createVariable('rain_stats','f8',('stats','number of timesteps'))

#########################################################

# mean 
phi_c_stats[0,:] = np.mean(np.mean(data_phi_c,axis=1),axis=0) # mean over ensemble members, and then over grid points
phi_stats[0,:] = np.mean(np.mean(data_phi,axis=1),axis=0) # mean over ensemble members, and then over grid points
wind_stats[0,:] = np.mean(np.mean(data_wind,axis=1),axis=0) # mean over ensemble members, and then over grid points
height_stats[0,:] = np.mean(np.mean(data_height,axis=1),axis=0) # mean over ensemble members, and then over grid points
rain_stats[0,:] = np.mean(np.mean(data_rain,axis=1),axis=0) # mean over ensemble members, and then over grid points

# median 
phi_c_stats[1,:] = np.mean(np.median(data_phi_c,axis=1),axis=0) # median over ensemble members, and then over grid points
phi_stats[1,:] = np.mean(np.median(data_phi,axis=1),axis=0) # median over ensemble members, and then over grid points
wind_stats[1,:] = np.mean(np.median(data_wind,axis=1),axis=0) # median over ensemble members, and then over grid points
height_stats[1,:] = np.mean(np.median(data_height,axis=1),axis=0) # median over ensemble members, and then over grid points
rain_stats[1,:] = np.mean(np.median(data_rain,axis=1),axis=0) # median over ensemble members, and then over grid points

# stdev 
phi_c_stats[2,:] = np.mean(np.std(data_phi_c,ddof=1,axis=1),axis=0) # stdev over ensemble members, and then over grid points
phi_stats[2,:] = np.mean(np.std(data_phi,ddof=1,axis=1),axis=0) # stdev over ensemble members, and then over grid points
wind_stats[2,:] = np.mean(np.std(data_wind,ddof=1,axis=1),axis=0) # stdev over ensemble members, and then over grid points
height_stats[2,:] = np.mean(np.std(data_height,ddof=1,axis=1),axis=0) # stdev over ensemble members, and then over grid points
rain_stats[2,:] = np.mean(np.std(data_rain,ddof=1,axis=1),axis=0) # stdev over ensemble members, and then over grid points

# absolute error 
phi_c_stats[3,:] = np.mean(np.mean(data_phi_c,axis=1) - data_phi_c_truth[:,0,:],axis=0)
phi_stats[3,:] = np.mean(np.mean(data_phi,axis=1) - data_phi_truth[:,0,:],axis=0)
wind_stats[3,:] = np.mean(np.mean(data_wind,axis=1) - data_wind_truth[:,0,:],axis=0)
height_stats[3,:] = np.mean(np.mean(data_height,axis=1) - data_height_truth[:,0,:],axis=0)
rain_stats[3,:] = np.mean(np.mean(data_rain,axis=1) - data_rain_truth[:,0,:],axis=0)

# normalised spread
phi_c_stats[4,:] = np.mean(np.std(data_phi_c,ddof=1,axis=1)/np.mean(data_phi_c,axis=1),axis=0) # stdev over ensemble members, and then over grid points
phi_stats[4,:] = np.mean(np.std(data_phi,ddof=1,axis=1)/np.mean(data_phi,axis=1),axis=0) # stdev over ensemble members, and then over grid points
wind_stats[4,:] = np.mean(np.std(data_wind,ddof=1,axis=1)/np.mean(data_wind,axis=1),axis=0) # stdev over ensemble members, and then over grid points
height_stats[4,:] = np.mean(np.std(data_height,ddof=1,axis=1)/np.mean(data_height,axis=1),axis=0) # stdev over ensemble members, and then over grid points
rain_stats[4,:] = np.mean(np.std(data_rain,ddof=1,axis=1)/np.mean(data_rain,axis=1),axis=0) # stdev over ensemble members, and then over grid points

# normalised spread 2
phi_c_stats[5,:] = phi_c_stats[2,:]/phi_c_stats[0,:]
phi_stats[5,:] = phi_stats[2,:]/phi_stats[0,:]
wind_stats[5,:] = wind_stats[2,:]/wind_stats[0,:]
height_stats[5,:] = height_stats[2,:]/height_stats[0,:]
rain_stats[5,:] = rain_stats[2,:]/rain_stats[0,:]

DIST.close()






