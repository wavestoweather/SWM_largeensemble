import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import random
import argparse
import numpy.ma as ma
from netCDF4 import Dataset
from datetime import datetime
import pickle

################################################

# creates convergence measure statistics for: 
# 1) mean
# 2) 95th percentile
# 3) standard deviation
# this is for the neighbourhood data!

# in this _2 version, it takes into the neighbourhood advantage properly. ie at an ensemble size of 10, it would have 10 * grid points in neighbourhood worth of ensemble members. 
# changed number and file saving names

################################################

parser = argparse.ArgumentParser(description='Outputs convergence measure. Input variables please.')
parser.add_argument('--grid_point',required=True, type=int, help='grid point')
parser.add_argument('--times',required=True, type=int, help='how many time points') # calculate for every hour after and including 5am
parser.add_argument('--times_diff',required=True, type=int, help='space between time points')
parser.add_argument('--times_start',required=True, type=int, help='start time for time points')
parser.add_argument('--boots',required=True, type=int, help='number of bootstrapped distributions to create')
parser.add_argument('--n_half_size',required=True, type=int, help='n_half_size') #input index of half size 
parser.add_argument('--step',required=True,type=int,help='')
args = parser.parse_args()
grid_point = args.grid_point
boots = args.boots
times = args.times
times_diff = args.times_diff
times_start = args.times_start
n_half_size = args.n_half_size
step = args.step 

#########################################################

# create function
def make_boot(data,boots):
    np.random.seed(24)
    bootstrap = np.zeros((boots,len(data)))
    for i in range(boots):
        print(i)
        bootstrap[i,:] = np.random.choice(data,size=len(data))
    return bootstrap
    
#########################################################

# import data
n_half_sizes = [2,10,20,50,100,200]

DIST = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/distribution_data_plus1/DIST_grid_point_500.nc','r')
wind_ = ma.getdata(DIST['/G_WIND/wind_dist_data'][n_half_size,:,:,:])
height_ = ma.getdata(DIST['/G_HEIGHT/height_dist_data'][n_half_size,:,:,:])
rain_ = ma.getdata(DIST['/G_RAIN/rain_dist_data'][n_half_size,:,:,:])

# prepare data
# wind 

wind_ = np.where(wind_<99999,wind_,np.nan)
height_ = np.where(height_<99999,height_,np.nan)
rain_ = np.where(rain_<99999,rain_,np.nan)

#########################################################

every = (n_half_sizes[n_half_size] * 2) + 1
every = every * step
number = 5000 * ((2 * n_half_sizes[n_half_size])+1) # number of max ensemble size you want to have in your convergence plot
points = np.int(number/every)
wind_mean_output = np.zeros((times,points))
height_mean_output = np.zeros((times,points))
rain_mean_output = np.zeros((times,points))
wind_q95_output = np.zeros((times,points))
height_q95_output = np.zeros((times,points))
rain_q95_output = np.zeros((times,points))
wind_stdev_output = np.zeros((times,points))
height_stdev_output = np.zeros((times,points))
rain_stdev_output = np.zeros((times,points))

# wind

for i in range(times):

 wind = wind_[:,:,times_start+(times_diff*i)][~ np.isnan(wind_[:,:,times_start+(times_diff*i)])]

 bootstrap = make_boot(wind,boots)

 # calculate for mean first
 mean_total = np.mean(wind)
 cumsum_boot = np.cumsum(bootstrap,axis=1)
 cumsum_boot_mean = np.zeros((boots,points))

 for j in range(boots):
     for k in range(points):
         cumsum_boot_mean[j,k] = cumsum_boot[j,(k+1)*every-1]/np.arange(every,number+1,every)[k]

 cumsum_boot_mean_rd = (cumsum_boot_mean-mean_total)

 # calculate width of CI
 wind_mean_output[i,:] = np.quantile(cumsum_boot_mean_rd,q=0.975,axis=0)-np.quantile(cumsum_boot_mean_rd,q=0.025,axis=0)

 # now calculate for 95th percentile!

 q95_total = np.quantile(wind,q=0.95)
 cumsum_boot_q95 = np.zeros((boots,points))

 for j in range(boots):
    for k in range(points):
        cumsum_boot_q95[j,k] = np.quantile(bootstrap[j,:every*(k+1)],q=0.95)

 cumsum_boot_q95_rd = cumsum_boot_q95 - q95_total

 # calculate width of CI
 wind_q95_output[i,:] = np.quantile(cumsum_boot_q95_rd,0.975,axis=0)-np.quantile(cumsum_boot_q95_rd,0.025,axis=0)

 # now calculate for standard deviation!

 stdev_total = np.std(wind)
 cumsum_boot_stdev = np.zeros((boots,points))

 for j in range(boots):
     for k in range(points):
         cumsum_boot_stdev[j,k] = np.std(bootstrap[j,:(k+1)*every])

 cumsum_boot_stdev_rd = cumsum_boot_stdev - stdev_total

  # calculate width of CI
 wind_stdev_output[i,:] = np.quantile(cumsum_boot_stdev_rd,0.975,axis=0)-np.quantile(cumsum_boot_stdev_rd,0.025,axis=0)

#save data
with open('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/convergence_data_plus1/convergence_data_neighbourhood_mean_weakforcing275000_gp_'+str(grid_point)+'half_size_'+str(n_half_sizes[n_half_size])+'_wind_time_start='+str(times_start)+'_2'+'step='+str(step),'wb') as f:
 pickle.dump(wind_mean_output,f)

with open('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/convergence_data_plus1/convergence_data_neighbourhood_q95_weakforcing275000_gp_'+str(grid_point)+'half_size_'+str(n_half_sizes[n_half_size])+'_wind_time_start='+str(times_start)+'_2'+'step='+str(step),'wb') as f:
 pickle.dump(wind_q95_output,f)

with open('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/convergence_data_plus1/convergence_data_neighbourhood_stdev_weakforcing275000_gp_'+str(grid_point)+'half_size_'+str(n_half_sizes[n_half_size])+'_wind_time_start='+str(times_start)+'_2'+'step='+str(step),'wb') as f:
 pickle.dump(wind_stdev_output,f)

#######################################################
# height

for i in range(times):

  height = height_[:,:,times_start+(times_diff*i)][~ np.isnan(height_[:,:,times_start+(times_diff*i)])]
  # create bootstrap

  bootstrap = make_boot(height,boots)

  # calculate for mean first
  mean_total = np.mean(height)
  cumsum_boot = np.cumsum(bootstrap,axis=1)
  cumsum_boot_mean = np.zeros((boots,points))

  for j in range(boots):
      for k in range(points):
          cumsum_boot_mean[j,k] = cumsum_boot[j,(k+1)*every-1]/np.arange(every,number+1,every)[k]

  cumsum_boot_mean_rd = (cumsum_boot_mean-mean_total)

#   calculate width of CI
  height_mean_output[i,:] = np.quantile(cumsum_boot_mean_rd,0.975,axis=0)-np.quantile(cumsum_boot_mean_rd,0.025,axis=0)

  # now calculate for 95th percentile!

  q95_total = np.quantile(height,q=0.95)
  cumsum_boot_q95 = np.zeros((boots,points))

  for j in range(boots):
      for k in range(points):
          cumsum_boot_q95[j,k] = np.quantile(bootstrap[j,:every*(k+1)],q=0.95)

  cumsum_boot_q95_rd = cumsum_boot_q95 - q95_total

 #   calculate width of CI
  height_q95_output[i,:] = np.quantile(cumsum_boot_q95_rd,0.975,axis=0)-np.quantile(cumsum_boot_q95_rd,0.025,axis=0)

#   now calculate for standard deviation!

  stdev_total = np.std(height)
  cumsum_boot_stdev = np.zeros((boots,points))

  for j in range(boots):
      for k in range(points):
          cumsum_boot_stdev[j,k] = np.std(bootstrap[j,:(k+1)*every])

  cumsum_boot_stdev_rd = cumsum_boot_stdev - stdev_total

#   calculate width of CI
  height_stdev_output[i,:] = np.quantile(cumsum_boot_stdev_rd,0.975,axis=0)-np.quantile(cumsum_boot_stdev_rd,0.025,axis=0)

# save data
with open('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/convergence_data_plus1/convergence_data_neighbourhood_mean_weakforcing275000_gp_'+str(grid_point)+'half_size_'+str(n_half_sizes[n_half_size])+'_height_time_start='+str(times_start)+'_2'+'step='+str(step),'wb') as f:
  pickle.dump(height_mean_output,f)

with open('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/convergence_data_plus1/convergence_data_neighbourhood_q95_weakforcing275000_gp_'+str(grid_point)+'half_size_'+str(n_half_sizes[n_half_size])+'_height_time_start='+str(times_start)+'_2'+'step='+str(step),'wb') as f:
  pickle.dump(height_q95_output,f)

with open('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/convergence_data_plus1/convergence_data_neighbourhood_stdev_weakforcing275000_gp_'+str(grid_point)+'half_size_'+str(n_half_sizes[n_half_size])+'_height_time_start='+str(times_start)+'_2'+'step='+str(step),'wb') as f:
  pickle.dump(height_stdev_output,f)


#########################################################
# rain

for i in range(times):

    rain = rain_[:,:,times_start+(times_diff*i)][~ np.isnan(rain_[:,:,times_start+(times_diff*i)])]
    # create bootstrap

    bootstrap = make_boot(rain,boots)

    # calculate for mean first
    mean_total = np.mean(rain)
    cumsum_boot = np.cumsum(bootstrap,axis=1)
    cumsum_boot_mean = np.zeros((boots,points))

    for j in range(boots):
        for k in range(points):
            cumsum_boot_mean[j,k] = cumsum_boot[j,(k+1)*every-1]/np.arange(every,number+1,every)[k]

    cumsum_boot_mean_rd = (cumsum_boot_mean-mean_total)

    # calculate width of CI
    rain_mean_output[i,:] = np.quantile(cumsum_boot_mean_rd,0.975,axis=0)-np.quantile(cumsum_boot_mean_rd,0.025,axis=0)

    # now calculate for 95th percentile!

    q95_total = np.quantile(rain,q=0.95)
    cumsum_boot_q95 = np.zeros((boots,points))

    for j in range(boots):
        for k in range(points):
            cumsum_boot_q95[j,k] = np.quantile(bootstrap[j,:every*(k+1)],q=0.95)

    cumsum_boot_q95_rd = cumsum_boot_q95 - q95_total

    # calculate width of CI
    rain_q95_output[i,:] = np.quantile(cumsum_boot_q95_rd,0.975,axis=0)-np.quantile(cumsum_boot_q95_rd,0.025,axis=0)

    # now calculate for standard deviation!

    stdev_total = np.std(rain)
    cumsum_boot_stdev = np.zeros((boots,points))

    for j in range(boots):
        for k in range(points):
            cumsum_boot_stdev[j,k] = np.std(bootstrap[j,:(k+1)*every])

    cumsum_boot_stdev_rd = cumsum_boot_stdev - stdev_total

    # calculate width of CI
    rain_stdev_output[i,:] = np.quantile(cumsum_boot_stdev_rd,0.975,axis=0)-np.quantile(cumsum_boot_stdev_rd,0.025,axis=0)

# save data
with open('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/convergence_data_plus1/convergence_data_neighbourhood_mean_weakforcing275000_gp_'+str(grid_point)+'half_size_'+str(n_half_sizes[n_half_size])+'_rain_time_start='+str(times_start)+'_2'+'step='+str(step),'wb') as f:
    pickle.dump(rain_mean_output,f)

with open('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/convergence_data_plus1/convergence_data_neighbourhood_q95_weakforcing275000_gp_'+str(grid_point)+'half_size_'+str(n_half_sizes[n_half_size])+'_rain_time_start='+str(times_start)+'_2'+'step='+str(step),'wb') as f:
    pickle.dump(rain_q95_output,f)

with open('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/convergence_data_plus1/convergence_data_neighbourhood_stdev_weakforcing275000_gp_'+str(grid_point)+'half_size_'+str(n_half_sizes[n_half_size])+'_rain_time_start='+str(times_start)+'_2'+'step='+str(step),'wb') as f:
    pickle.dump(rain_stdev_output,f)
