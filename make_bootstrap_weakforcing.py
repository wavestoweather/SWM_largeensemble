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

################################################

parser = argparse.ArgumentParser(description='Outputs convergence measure. Input variables please.')
parser.add_argument('--grid_point',required=True, type=int, help='grid point')
parser.add_argument('--times',required=True, type=int, help='how many time points') # calculate for every hour after and including 5am
parser.add_argument('--times_diff',required=True, type=int, help='space between time points')
parser.add_argument('--times_start',required=True, type=int, help='start time for time points')
parser.add_argument('--boots',required=True, type=int, help='number of bootstrapped distributions to create')
parser.add_argument('--every',required=True, type=int, help='sample at every ensemble size possible?')
args = parser.parse_args()
grid_point = args.grid_point
boots = args.boots
every = args.every
times = args.times
times_diff = args.times_diff
times_start = args.times_start

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
wind = pickle.load(open('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/distribution_data/distribution_data_weakforcing275000_gp_'+str(grid_point)+'_wind','rb'))
height = pickle.load(open('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/distribution_data/distribution_data_weakforcing275000_gp_'+str(grid_point)+'_height','rb'))
rain = pickle.load(open('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/distribution_data/distribution_data_weakforcing275000_gp_'+str(grid_point)+'_rain','rb'))

#########################################################

number = 5000 # number of max ensemble size you want
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
    # create bootstrap

    bootstrap = make_boot(wind[:,times_start+(times_diff*i)],boots)

    # calculate for mean first
    mean_total = np.mean(wind[:,times_start+(times_diff*i)])
    cumsum_boot = np.cumsum(bootstrap,axis=1)
    cumsum_boot_mean = np.zeros((boots,points))

    for j in range(boots):
        for k in range(points):
            cumsum_boot_mean[j,k] = cumsum_boot[j,(k+1)*every-1]/np.arange(every,number+1,every)[k]

    cumsum_boot_mean_rd = (cumsum_boot_mean-mean_total)

    # calculate width of CI
    wind_mean_output[i,:] = np.quantile(cumsum_boot_mean_rd,q=0.975,axis=0)-np.quantile(cumsum_boot_mean_rd,q=0.025,axis=0)

    # now calculate for 95th percentile!

    q95_total = np.quantile(wind[:,times_start+(times_diff*i)],q=0.95)
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

# save data
with open('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/convergence_data/convergence_data_mean_weakforcing275000_gp_'+str(grid_point)+'times_start='+str(times_start)+'_wind_every='+str(every),'wb') as f:
    pickle.dump(wind_mean_output,f)

with open('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/convergence_data/convergence_data_q95_weakforcing275000_gp_'+str(grid_point)+'times_start='+str(times_start)+'_wind_every='+str(every),'wb') as f:
    pickle.dump(wind_q95_output,f)

with open('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/convergence_data/convergence_data_stdev_weakforcing275000_gp_'+str(grid_point)+'times_start='+str(times_start)+'_wind_every='+str(every),'wb') as f:
    pickle.dump(wind_stdev_output,f)

#########################################################
# height

for i in range(times):
    # create bootstrap

    bootstrap = make_boot(height[:,times_start+(times_diff*i)],boots)

    # calculate for mean first
    mean_total = np.mean(height[:,times_start+(times_diff*i)])
    cumsum_boot = np.cumsum(bootstrap,axis=1)
    cumsum_boot_mean = np.zeros((boots,points))

    for j in range(boots):
        for k in range(points):
            cumsum_boot_mean[j,k] = cumsum_boot[j,(k+1)*every-1]/np.arange(every,number+1,every)[k]

    cumsum_boot_mean_rd = (cumsum_boot_mean-mean_total)

    # calculate width of CI
    height_mean_output[i,:] = np.quantile(cumsum_boot_mean_rd,0.975,axis=0)-np.quantile(cumsum_boot_mean_rd,0.025,axis=0)

    # now calculate for 95th percentile!

    q95_total = np.quantile(height[:,times_start+(times_diff*i)],q=0.95)
    cumsum_boot_q95 = np.zeros((boots,points))

    for j in range(boots):
        for k in range(points):
            cumsum_boot_q95[j,k] = np.quantile(bootstrap[j,:every*(k+1)],q=0.95)

    cumsum_boot_q95_rd = cumsum_boot_q95 - q95_total

    # calculate width of CI
    height_q95_output[i,:] = np.quantile(cumsum_boot_q95_rd,0.975,axis=0)-np.quantile(cumsum_boot_q95_rd,0.025,axis=0)

    # now calculate for standard deviation!

    stdev_total = np.std(height)
    cumsum_boot_stdev = np.zeros((boots,points))

    for j in range(boots):
        for k in range(points):
            cumsum_boot_stdev[j,k] = np.std(bootstrap[j,:(k+1)*every])

    cumsum_boot_stdev_rd = cumsum_boot_stdev - stdev_total

    # calculate width of CI
    height_stdev_output[i,:] = np.quantile(cumsum_boot_stdev_rd,0.975,axis=0)-np.quantile(cumsum_boot_stdev_rd,0.025,axis=0)

# save data
with open('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/convergence_data/convergence_data_mean_weakforcing275000_gp_'+str(grid_point)+'times_start='+str(times_start)+'_height_every='+str(every),'wb') as f:
    pickle.dump(height_mean_output,f)

with open('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/convergence_data/convergence_data_q95_weakforcing275000_gp_'+str(grid_point)+'times_start='+str(times_start)+'_height_every='+str(every),'wb') as f:
    pickle.dump(height_q95_output,f)

with open('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/convergence_data/convergence_data_stdev_weakforcing275000_gp_'+str(grid_point)+'times_start='+str(times_start)+'_height_every='+str(every),'wb') as f:
    pickle.dump(height_stdev_output,f)


#########################################################
# rain

for i in range(times):
    # create bootstrap

    bootstrap = make_boot(rain[:,times_start+(times_diff*i)],boots)

    # calculate for mean first
    mean_total = np.mean(rain[:,times_start+(times_diff*i)])
    cumsum_boot = np.cumsum(bootstrap,axis=1)
    cumsum_boot_mean = np.zeros((boots,points))

    for j in range(boots):
        for k in range(points):
            cumsum_boot_mean[j,k] = cumsum_boot[j,(k+1)*every-1]/np.arange(every,number+1,every)[k]

    cumsum_boot_mean_rd = (cumsum_boot_mean-mean_total)

    # calculate width of CI
    rain_mean_output[i,:] = np.quantile(cumsum_boot_mean_rd,0.975,axis=0)-np.quantile(cumsum_boot_mean_rd,0.025,axis=0)

    # now calculate for 95th percentile!

    q95_total = np.quantile(rain[:,times_start+(times_diff*i)],q=0.95)
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
with open('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/convergence_data/convergence_data_mean_weakforcing275000_gp_'+str(grid_point)+'times_start='+str(times_start)+'_rain_every='+str(every),'wb') as f:
    pickle.dump(rain_mean_output,f)

with open('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/convergence_data/convergence_data_q95_weakforcing275000_gp_'+str(grid_point)+'times_start='+str(times_start)+'_rain_every='+str(every),'wb') as f:
    pickle.dump(rain_q95_output,f)

with open('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/convergence_data/convergence_data_stdev_weakforcing275000_gp_'+str(grid_point)+'times_start='+str(times_start)+'_rain_every='+str(every),'wb') as f:
    pickle.dump(rain_stdev_output,f)
