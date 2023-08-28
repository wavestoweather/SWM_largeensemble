import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from netCDF4 import Dataset
import numpy.ma as ma
from datetime import datetime
import seaborn as sns
import matplotlib
import os
import argparse
import scipy.stats as ss
import pickle

sns.set(color_codes = True)
sns.set_style('whitegrid')

################################################

# plot used to compare convergence between the strong and weak forcings. Store in weak_forcing_1's directories.

# change when using script: 
# file for input data

################################################

parser = argparse.ArgumentParser(description='Plots comparison convergence plots with an emphasis on estimating the effective ensemble size (ees). Input variables please.')
parser.add_argument('--time',required=True, type=int, help='time') 
args = parser.parse_args() 
time = args.time

#################################################################################################

# load data

n_half_sizes = [2,10,20,50,100,200,500]

# weak forcing 

wind_mean_2 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[0])+'_wind_2','rb'))
wind_stdev_2 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[0])+'_wind_2','rb'))
wind_q95_2 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[0])+'_wind_2','rb'))
wind_mean_10 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[1])+'_wind_2','rb'))
wind_stdev_10 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[1])+'_wind_2','rb'))
wind_q95_10 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[1])+'_wind_2','rb'))
wind_mean_20 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[2])+'_wind_2','rb'))
wind_stdev_20 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[2])+'_wind_2','rb'))
wind_q95_20 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[2])+'_wind_2','rb'))
wind_mean_50 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[3])+'_wind_2','rb'))
wind_stdev_50 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[3])+'_wind_2','rb'))
wind_q95_50 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[3])+'_wind_2','rb'))
wind_mean_100 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[4])+'_wind_2','rb'))
wind_stdev_100 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[4])+'_wind_2','rb'))
wind_q95_100 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[4])+'_wind_2','rb'))
wind_mean_200 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[5])+'_wind_2','rb'))
wind_stdev_200 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[5])+'_wind_2','rb'))
wind_q95_200 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[5])+'_wind_2','rb'))
#wind_mean_500 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[6])+'_wind_2','rb'))
#wind_stdev_500 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[6])+'_wind_2','rb'))
#wind_q95_500 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[6])+'_wind_2','rb'))

height_mean_2 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[0])+'_height_2','rb'))
height_stdev_2 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[0])+'_height_2','rb'))
height_q95_2 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[0])+'_height_2','rb'))
height_mean_10 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[1])+'_height_2','rb'))
height_stdev_10 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[1])+'_height_2','rb'))
height_q95_10 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[1])+'_height_2','rb'))
height_mean_20 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[2])+'_height_2','rb'))
height_stdev_20 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[2])+'_height_2','rb'))
height_q95_20 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[2])+'_height_2','rb'))
height_mean_50 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[3])+'_height_2','rb'))
height_stdev_50 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[3])+'_height_2','rb'))
height_q95_50 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[3])+'_height_2','rb'))
height_mean_100 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[4])+'_height_2','rb'))
height_stdev_100 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[4])+'_height_2','rb'))
height_q95_100 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[4])+'_height_2','rb'))
height_mean_200 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[5])+'_height_2','rb'))
height_stdev_200 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[5])+'_height_2','rb'))
height_q95_200 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[5])+'_height_2','rb'))
#height_mean_500 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[6])+'_height_2','rb'))
#height_stdev_500 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[6])+'_height_2','rb'))
#height_q95_500 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[6])+'_height_2','rb'))

rain_mean_2 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[0])+'_rain_2','rb'))
rain_stdev_2 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[0])+'_rain_2','rb'))
rain_q95_2 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[0])+'_rain_2','rb'))
rain_mean_10 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[1])+'_rain_2','rb'))
rain_stdev_10 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[1])+'_rain_2','rb'))
rain_q95_10 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[1])+'_rain_2','rb'))
rain_mean_20 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[2])+'_rain_2','rb'))
rain_stdev_20 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[2])+'_rain_2','rb'))
rain_q95_20 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[2])+'_rain_2','rb'))
rain_mean_50 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[3])+'_rain_2','rb'))
rain_stdev_50 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[3])+'_rain_2','rb'))
rain_q95_50 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[3])+'_rain_2','rb'))
rain_mean_100 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[4])+'_rain_2','rb'))
rain_stdev_100 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[4])+'_rain_2','rb'))
rain_q95_100 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[4])+'_rain_2','rb'))
rain_mean_200 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[5])+'_rain_2','rb'))
rain_stdev_200 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[5])+'_rain_2','rb'))
rain_q95_200 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[5])+'_rain_2','rb'))
#rain_mean_500 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[6])+'_rain_2','rb'))
#rain_stdev_500 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[6])+'_rain_2','rb'))
#rain_q95_500 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_weakforcing_gp_500half_size_'+str(n_half_sizes[6])+'_rain_2','rb'))

# import single grid point data 
wind_mean = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_mean_extendedrun1_weakforcing_gp_500_wind','rb'))
height_mean = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_mean_extendedrun1_weakforcing_gp_500_height','rb'))
rain_mean = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_mean_extendedrun1_weakforcing_gp_500_rain','rb'))
wind_stdev = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_stdev_extendedrun1_weakforcing_gp_500_wind','rb'))
height_stdev = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_stdev_extendedrun1_weakforcing_gp_500_height','rb'))
rain_stdev = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_stdev_extendedrun1_weakforcing_gp_500_rain','rb'))
wind_q95 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_q95_extendedrun1_weakforcing_gp_500_wind','rb'))
height_q95 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_q95_extendedrun1_weakforcing_gp_500_height','rb'))
rain_q95 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_data_q95_extendedrun1_weakforcing_gp_500_rain','rb'))

# strong forcing 

sf_wind_mean_2 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[0])+'_wind_2','rb'))
sf_wind_stdev_2 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[0])+'_wind_2','rb'))
sf_wind_q95_2 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[0])+'_wind_2','rb'))
sf_wind_mean_10 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[1])+'_wind_2','rb'))
sf_wind_stdev_10 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[1])+'_wind_2','rb'))
sf_wind_q95_10 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[1])+'_wind_2','rb'))
sf_wind_mean_20 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[2])+'_wind_2','rb'))
sf_wind_stdev_20 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[2])+'_wind_2','rb'))
sf_wind_q95_20 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[2])+'_wind_2','rb'))
sf_wind_mean_50 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[3])+'_wind_2','rb'))
sf_wind_stdev_50 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[3])+'_wind_2','rb'))
sf_wind_q95_50 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[3])+'_wind_2','rb'))
sf_wind_mean_100 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[4])+'_wind_2','rb'))
sf_wind_stdev_100 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[4])+'_wind_2','rb'))
sf_wind_q95_100 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[4])+'_wind_2','rb'))
sf_wind_mean_200 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[5])+'_wind_2','rb'))
sf_wind_stdev_200 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[5])+'_wind_2','rb'))
sf_wind_q95_200 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[5])+'_wind_2','rb'))
#sf_wind_mean_500 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[6])+'_wind_2','rb'))
#sf_wind_stdev_500 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[6])+'_wind_2','rb'))
#sf_wind_q95_500 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[6])+'_wind_2','rb'))

sf_height_mean_2 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[0])+'_height_2','rb'))
sf_height_stdev_2 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[0])+'_height_2','rb'))
sf_height_q95_2 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[0])+'_height_2','rb'))
sf_height_mean_10 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[1])+'_height_2','rb'))
sf_height_stdev_10 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[1])+'_height_2','rb'))
sf_height_q95_10 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[1])+'_height_2','rb'))
sf_height_mean_20 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[2])+'_height_2','rb'))
sf_height_stdev_20 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[2])+'_height_2','rb'))
sf_height_q95_20 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[2])+'_height_2','rb'))
sf_height_mean_50 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[3])+'_height_2','rb'))
sf_height_stdev_50 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[3])+'_height_2','rb'))
sf_height_q95_50 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[3])+'_height_2','rb'))
sf_height_mean_100 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[4])+'_height_2','rb'))
sf_height_stdev_100 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[4])+'_height_2','rb'))
sf_height_q95_100 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[4])+'_height_2','rb'))
sf_height_mean_200 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[5])+'_height_2','rb'))
sf_height_stdev_200 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[5])+'_height_2','rb'))
sf_height_q95_200 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[5])+'_height_2','rb'))
#sf_height_mean_500 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[6])+'_height_2','rb'))
#sf_height_stdev_500 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[6])+'_height_2','rb'))
#sf_height_q95_500 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[6])+'_height_2','rb'))

sf_rain_mean_2 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[0])+'_rain_2','rb'))
sf_rain_stdev_2 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[0])+'_rain_2','rb'))
sf_rain_q95_2 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[0])+'_rain_2','rb'))
sf_rain_mean_10 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[1])+'_rain_2','rb'))
sf_rain_stdev_10 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[1])+'_rain_2','rb'))
sf_rain_q95_10 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[1])+'_rain_2','rb'))
sf_rain_mean_20 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[2])+'_rain_2','rb'))
sf_rain_stdev_20 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[2])+'_rain_2','rb'))
sf_rain_q95_20 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[2])+'_rain_2','rb'))
sf_rain_mean_50 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[3])+'_rain_2','rb'))
sf_rain_stdev_50 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[3])+'_rain_2','rb'))
sf_rain_q95_50 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[3])+'_rain_2','rb'))
sf_rain_mean_100 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[4])+'_rain_2','rb'))
sf_rain_stdev_100 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[4])+'_rain_2','rb'))
sf_rain_q95_100 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[4])+'_rain_2','rb'))
sf_rain_mean_200 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[5])+'_rain_2','rb'))
sf_rain_stdev_200 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[5])+'_rain_2','rb'))
sf_rain_q95_200 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[5])+'_rain_2','rb'))
#sf_rain_mean_500 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_mean_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[6])+'_rain_2','rb'))
#sf_rain_stdev_500 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_stdev_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[6])+'_rain_2','rb'))
#sf_rain_q95_500 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_neighbourhood_q95_extendedrun1_strongforcing15_gp_500half_size_'+str(n_half_sizes[6])+'_rain_2','rb'))

# import single grid point data 
sf_wind_mean = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_mean_extendedrun1_strongforcing15_singlegp_500_wind','rb'))
sf_height_mean = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_mean_extendedrun1_strongforcing15_singlegp_500_height','rb'))
sf_rain_mean = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_mean_extendedrun1_strongforcing15_singlegp_500_rain','rb'))
sf_wind_stdev = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_stdev_extendedrun1_strongforcing15_singlegp_500_wind','rb'))
sf_height_stdev = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_stdev_extendedrun1_strongforcing15_singlegp_500_height','rb'))
sf_rain_stdev = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_stdev_extendedrun1_strongforcing15_singlegp_500_rain','rb'))
sf_wind_q95 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_q95_extendedrun1_strongforcing15_singlegp_500_wind','rb'))
sf_height_q95 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_q95_extendedrun1_strongforcing15_singlegp_500_height','rb'))
sf_rain_q95 = pickle.load(open('/scratch/k/K.Tempest/Model_extension/strong_forcing_15/convergence_data/convergence_data_q95_extendedrun1_strongforcing15_singlegp_500_rain','rb'))

# define function and calculate b variable

x = np.arange(1,121)

def func(x,b):
	y = (10**b)*x**(-1/2)
	return y 

wind_mean_b = np.sum(((np.log10(wind_mean_20))) + (0.5 * np.log10(x)),axis=1)/len(wind_mean_20[0,:]) # change this to 200 when have full data!!!
height_mean_b = np.sum(((np.log10(height_mean_20))) + (0.5 * np.log10(x)),axis=1)/len(height_mean_20[0,:])
rain_mean_b = np.sum(((np.log10(rain_mean_20))) + (0.5 * np.log10(x)),axis=1)/len(rain_mean_20[0,:])

wind_q95_b = np.sum(((np.log10(wind_q95_20))) + (0.5 * np.log10(x)),axis=1)/len(wind_q95_20[0,:])
height_q95_b = np.sum(((np.log10(height_q95_20))) + (0.5 * np.log10(x)),axis=1)/len(height_q95_20[0,:])
rain_q95_b = np.sum(((np.log10(rain_q95_20))) + (0.5 * np.log10(x)),axis=1)/len(rain_q95_20[0,:])

wind_stdev_b = np.sum(((np.log10(wind_stdev_20[:,1:]))) + (0.5 * np.log10(x[1:])),axis=1)/len(wind_stdev_20[0,1:])
height_stdev_b = np.sum(((np.log10(height_stdev_20[:,1:]))) + (0.5 * np.log10(x[1:])),axis=1)/len(height_stdev_20[0,1:])
rain_stdev_b = np.sum(((np.log10(rain_stdev_20[:,1:]))) + (0.5 * np.log10(x[1:])),axis=1)/len(rain_stdev_20[0,1:])

# plot mean #

fig, (ax1, ax2, ax3)= plt.subplots(1, 3,figsize=[14,5])
plt.subplots_adjust(wspace = 0.5,hspace = 0)
plt.grid(True)
#################################################################################################

# do for only one time

ax1.plot(x,sf_wind_mean[time],linestyle='-',lw=4,alpha=0.5,color = 'dodgerblue',label='thicker=strong forcing 15') # strong forcing 
ax1.plot(x,sf_wind_mean_2[time],linestyle='-',lw=4,alpha=0.5,color = 'firebrick')
ax1.plot(x,sf_wind_mean_10[time],linestyle='-',lw=4,alpha=0.5,color = 'darkorange')
ax1.plot(x,sf_wind_mean_20[time],linestyle='-',lw=4,alpha=0.5,color = 'goldenrod')
ax1.plot(x,sf_wind_mean_50[time],linestyle='-',lw=4,alpha=0.5,color = 'chartreuse')
ax1.plot(x,sf_wind_mean_100[time],linestyle='-',lw=4,alpha=0.5,color = 'forestgreen')
ax1.plot(x,sf_wind_mean_200[time],linestyle='-',lw=4,alpha=0.5,color = 'rosybrown') 
#ax1.plot(x,sf_wind_mean_500[time],linestyle='-',lw=4,alpha=0.5,color = 'mediumturquoise')  

f = func(x, wind_mean_b[time]) #weak forcing 
ax1.fill_between(x,f+0.05*f,f-0.05*f,color='silver',alpha=0.5)
ax1.plot(x,wind_mean[time],linestyle='-',lw=1,color = 'dodgerblue',label='single gp')
ax1.plot(x,wind_mean_2[time],linestyle='-',lw=1,color = 'firebrick',label='r=2')
ax1.plot(x,wind_mean_10[time],linestyle='-',lw=1,color = 'darkorange',label='r=10')
ax1.plot(x,wind_mean_20[time],linestyle='-',lw=1,color = 'goldenrod',label='r=20')
ax1.plot(x,wind_mean_50[time],linestyle='-',lw=1,color = 'chartreuse',label='r=50')
ax1.plot(x,wind_mean_100[time],linestyle='-',lw=1,color = 'forestgreen',label='r=100')
ax1.plot(x,wind_mean_200[time],linestyle='-',lw=1,color = 'rosybrown',label='r=200 i='+str(time+5)+'hours a='+str("{:.2f}".format(10**wind_mean_b[time]))) 
#ax1.plot(x,wind_mean_500[time],linestyle='-',lw=1,color = 'mediumturquoise',label='r=500')  
    
ax1.legend()
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.tick_params(axis='x', labelsize=19)
ax1.tick_params(axis='y', labelsize=19)
ax1.set_ylabel('$Width$ $of$ $95\%$ $CI$ $(ms^{-1})$',fontsize=19);
ax1.set_xlabel('$Ensemble$ $size$ $(n)$',fontsize=19);
ax1.set_title('$Wind$ $mean$',fontsize=19)

#################################################################################################
ax2.plot(x,sf_height_mean[time],linestyle='-',lw=4,alpha=0.5,color = 'dodgerblue',label='thicker=strong forcing 15') # strong forcing 
ax2.plot(x,sf_height_mean_2[time],linestyle='-',lw=4,alpha=0.5,color = 'firebrick')
ax2.plot(x,sf_height_mean_10[time],linestyle='-',lw=4,alpha=0.5,color = 'darkorange')
ax2.plot(x,sf_height_mean_20[time],linestyle='-',lw=4,alpha=0.5,color = 'goldenrod')
ax2.plot(x,sf_height_mean_50[time],linestyle='-',lw=4,alpha=0.5,color = 'chartreuse')
ax2.plot(x,sf_height_mean_100[time],linestyle='-',lw=4,alpha=0.5,color = 'forestgreen')
ax2.plot(x,sf_height_mean_200[time],linestyle='-',lw=4,alpha=0.5,color = 'rosybrown') 
#ax2.plot(x,sf_height_mean_500[time],linestyle='-',lw=4,alpha=0.5,color = 'mediumturquoise') 

f = func(x, height_mean_b[time])
ax2.fill_between(x,f+0.05*f,f-0.05*f,color='silver',alpha=0.5)
ax2.plot(x,height_mean[time],linestyle='-',lw=1,color = 'dodgerblue',label='single gp')
ax2.plot(x,height_mean_2[time],linestyle='-',lw=1,color = 'firebrick',label='r=2')
ax2.plot(x,height_mean_10[time],linestyle='-',lw=1,color = 'darkorange',label='r=10')
ax2.plot(x,height_mean_20[time],linestyle='-',lw=1,color = 'goldenrod',label='r=20')
ax2.plot(x,height_mean_50[time],linestyle='-',lw=1,color = 'chartreuse',label='r=50')
ax2.plot(x,height_mean_100[time],linestyle='-',lw=1,color = 'forestgreen',label='r=100')
ax2.plot(x,height_mean_200[time],linestyle='-',lw=1,color = 'rosybrown',label='r=200 i='+str(time+5)+'hours a='+str("{:.2f}".format(10**height_mean_b[time]))) 
#ax2.plot(x,height_mean_500[time],linestyle='-',lw=1,color = 'mediumturquoise',label='r=500')  
    
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.tick_params(axis='x', labelsize=19)
ax2.tick_params(axis='y', labelsize=19)
ax2.set_ylabel('$Width$ $of$ $95\%$ $CI$ $(m)$',fontsize=19);
ax2.set_xlabel('$Ensemble$ $size$ $(n)$',fontsize=19);
ax2.set_title('$Height$ $mean$',fontsize=19)

#################################################################################################

ax3.plot(x,sf_rain_mean[time],linestyle='-',lw=4,alpha=0.5,color = 'dodgerblue',label='thicker=strong forcing 15') # strong forcing 
ax3.plot(x,sf_rain_mean_2[time],linestyle='-',lw=4,alpha=0.5,color = 'firebrick')
ax3.plot(x,sf_rain_mean_10[time],linestyle='-',lw=4,alpha=0.5,color = 'darkorange')
ax3.plot(x,sf_rain_mean_20[time],linestyle='-',lw=4,alpha=0.5,color = 'goldenrod')
ax3.plot(x,sf_rain_mean_50[time],linestyle='-',lw=4,alpha=0.5,color = 'chartreuse')
ax3.plot(x,sf_rain_mean_100[time],linestyle='-',lw=4,alpha=0.5,color = 'forestgreen')
ax3.plot(x,sf_rain_mean_200[time],linestyle='-',lw=4,alpha=0.5,color = 'rosybrown') 
#ax3.plot(x,sf_rain_mean_500[time],linestyle='-',lw=4,alpha=0.5,color = 'mediumturquoise') 

f = func(x, rain_mean_b[time])
ax3.fill_between(x,f+0.05*f,f-0.05*f,color='silver',alpha=0.5)
ax3.plot(x,rain_mean[time],linestyle='-',lw=1,color = 'dodgerblue',label='single gp')
ax3.plot(x,rain_mean_2[time],linestyle='-',lw=1,color = 'firebrick',label='r=2')
ax3.plot(x,rain_mean_10[time],linestyle='-',lw=1,color = 'darkorange',label='r=10')
ax3.plot(x,rain_mean_20[time],linestyle='-',lw=1,color = 'goldenrod',label='r=20')
ax3.plot(x,rain_mean_50[time],linestyle='-',lw=1,color = 'chartreuse',label='r=50')
ax3.plot(x,rain_mean_100[time],linestyle='-',lw=1,color = 'firebrick',label='r=100')
ax3.plot(x,rain_mean_200[time],linestyle='-',lw=1,color = 'forestgreen',label='r=200 i='+str(time+5)+'hours a='+str("{:.2f}".format(10**rain_mean_b[time]))) 
#ax3.plot(x,rain_mean_500[time],linestyle='-',lw=1,color = 'mediumturquoise',label='r=500')  
    
ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.tick_params(axis='x', labelsize=19)
ax3.tick_params(axis='y', labelsize=19)
ax3.set_ylabel('$Width$ $of$ $95\%$ $CI$',fontsize=19);
ax3.set_xlabel('$Ensemble$ $size$ $(n)$',fontsize=19);
ax3.set_title('$Rain$ $mean$',fontsize=19)

#################################################################################################
plt.tight_layout()
plt.savefig('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_plot_neighbourhood_eescomparison_mean_extendedrun1_weakforcing_time_'+str(time)+'.pdf')

# plot stdev #

fig, (ax1, ax2, ax3)= plt.subplots(1, 3,figsize=[14,5])
plt.subplots_adjust(wspace = 0.5,hspace = 0)
plt.grid(True)
#################################################################################################

# do for only one time

ax1.plot(x,sf_wind_stdev[time],linestyle='-',lw=4,alpha=0.5,color = 'dodgerblue',label='thicker=strong forcing 15') # strong forcing 
ax1.plot(x,sf_wind_stdev_2[time],linestyle='-',lw=4,alpha=0.5,color = 'firebrick')
ax1.plot(x,sf_wind_stdev_10[time],linestyle='-',lw=4,alpha=0.5,color = 'darkorange')
ax1.plot(x,sf_wind_stdev_20[time],linestyle='-',lw=4,alpha=0.5,color = 'goldenrod')
ax1.plot(x,sf_wind_stdev_50[time],linestyle='-',lw=4,alpha=0.5,color = 'chartreuse')
ax1.plot(x,sf_wind_stdev_100[time],linestyle='-',lw=4,alpha=0.5,color = 'forestgreen')
ax1.plot(x,sf_wind_stdev_200[time],linestyle='-',lw=4,alpha=0.5,color = 'rosybrown') 
#ax1.plot(x,sf_wind_stdev_500[time],linestyle='-',lw=4,alpha=0.5,color = 'mediumturquoise')  

f = func(x, wind_stdev_b[time])
ax1.fill_between(x,f+0.05*f,f-0.05*f,color='silver',alpha=0.5)
ax1.plot(x,wind_stdev[time],linestyle='-',lw=1,color = 'dodgerblue',label='single gp')
ax1.plot(x,wind_stdev_2[time],linestyle='-',lw=1,color = 'firebrick',label='r=2')
ax1.plot(x,wind_stdev_10[time],linestyle='-',lw=1,color = 'darkorange',label='r=10')
ax1.plot(x,wind_stdev_20[time],linestyle='-',lw=1,color = 'goldenrod',label='r=20')
ax1.plot(x,wind_stdev_50[time],linestyle='-',lw=1,color = 'chartreuse',label='r=50')
ax1.plot(x,wind_stdev_100[time],linestyle='-',lw=1,color = 'forestgreen',label='r=100')
ax1.plot(x,wind_stdev_200[time],linestyle='-',lw=1,color = 'rosybrown',label='r=200 i='+str(time+5)+'hours a='+str("{:.2f}".format(10**wind_stdev_b[time]))) 
#ax1.plot(x,wind_stdev_500[time],linestyle='-',lw=1,color = 'mediumturquoise',label='r=500')  
    
ax1.legend()
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.tick_params(axis='x', labelsize=19)
ax1.tick_params(axis='y', labelsize=19)
ax1.set_ylabel('$Width$ $of$ $95\%$ $CI$',fontsize=19);
ax1.set_xlabel('$Ensemble$ $size$ $(n)$',fontsize=19);
ax1.set_title('$Wind$ $stdev$',fontsize=19)

#################################################################################################

ax2.plot(x,sf_height_stdev[time],linestyle='-',lw=4,alpha=0.5,color = 'dodgerblue',label='thicker=strong forcing 15') # strong forcing 
ax2.plot(x,sf_height_stdev_2[time],linestyle='-',lw=4,alpha=0.5,color = 'firebrick')
ax2.plot(x,sf_height_stdev_10[time],linestyle='-',lw=4,alpha=0.5,color = 'darkorange')
ax2.plot(x,sf_height_stdev_20[time],linestyle='-',lw=4,alpha=0.5,color = 'goldenrod')
ax2.plot(x,sf_height_stdev_50[time],linestyle='-',lw=4,alpha=0.5,color = 'chartreuse')
ax2.plot(x,sf_height_stdev_100[time],linestyle='-',lw=4,alpha=0.5,color = 'forestgreen')
ax2.plot(x,sf_height_stdev_200[time],linestyle='-',lw=4,alpha=0.5,color = 'rosybrown') 
#ax2.plot(x,sf_height_stdev_500[time],linestyle='-',lw=4,alpha=0.5,color = 'mediumturquoise')  

f = func(x, height_stdev_b[time])
ax2.fill_between(x,f+0.05*f,f-0.05*f,color='silver',alpha=0.5)
ax2.plot(x,height_stdev[time],linestyle='-',lw=1,color = 'dodgerblue',label='single gp')
ax2.plot(x,height_stdev_2[time],linestyle='-',lw=1,color = 'firebrick',label='r=2')
ax2.plot(x,height_stdev_10[time],linestyle='-',lw=1,color = 'darkorange',label='r=10')
ax2.plot(x,height_stdev_20[time],linestyle='-',lw=1,color = 'goldenrod',label='r=20')
ax2.plot(x,height_stdev_50[time],linestyle='-',lw=1,color = 'chartreuse',label='r=50')
ax2.plot(x,height_stdev_100[time],linestyle='-',lw=1,color = 'forestgreen',label='r=100')
ax2.plot(x,height_stdev_200[time],linestyle='-',lw=1,color = 'rosybrown',label='r=200 i='+str(time+5)+'hours a='+str("{:.2f}".format(10**height_stdev_b[time]))) 
#ax2.plot(x,height_stdev_500[time],linestyle='-',lw=1,color = 'mediumturquoise',label='r=500')  
    
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.tick_params(axis='x', labelsize=19)
ax2.tick_params(axis='y', labelsize=19)
ax2.set_ylabel('$Width$ $of$ $95\%$ $CI$',fontsize=19);
ax2.set_xlabel('$Ensemble$ $size$ $(n)$',fontsize=19);
ax2.set_title('$Height$ $stdev$',fontsize=19)

#################################################################################################

ax3.plot(x,sf_rain_stdev[time],linestyle='-',lw=4,alpha=0.5,color = 'dodgerblue',label='thicker=strong forcing 15') # strong forcing 
ax3.plot(x,sf_rain_stdev_2[time],linestyle='-',lw=4,alpha=0.5,color = 'firebrick')
ax3.plot(x,sf_rain_stdev_10[time],linestyle='-',lw=4,alpha=0.5,color = 'darkorange')
ax3.plot(x,sf_rain_stdev_20[time],linestyle='-',lw=4,alpha=0.5,color = 'goldenrod')
ax3.plot(x,sf_rain_stdev_50[time],linestyle='-',lw=4,alpha=0.5,color = 'chartreuse')
ax3.plot(x,sf_rain_stdev_100[time],linestyle='-',lw=4,alpha=0.5,color = 'forestgreen')
ax3.plot(x,sf_rain_stdev_200[time],linestyle='-',lw=4,alpha=0.5,color = 'rosybrown') 
#ax3.plot(x,sf_rain_stdev_500[time],linestyle='-',lw=4,alpha=0.5,color = 'mediumturquoise')  

f = func(x, rain_stdev_b[time])
ax3.fill_between(x,f+0.05*f,f-0.05*f,color='silver',alpha=0.5)
ax3.plot(x,rain_stdev[time],linestyle='-',lw=1,color = 'dodgerblue',label='single gp')
ax3.plot(x,rain_stdev_2[time],linestyle='-',lw=1,color = 'firebrick',label='r=2')
ax3.plot(x,rain_stdev_10[time],linestyle='-',lw=1,color = 'darkorange',label='r=10')
ax3.plot(x,rain_stdev_20[time],linestyle='-',lw=1,color = 'goldenrod',label='r=20')
ax3.plot(x,rain_stdev_50[time],linestyle='-',lw=1,color = 'chartreuse',label='r=50')
ax3.plot(x,rain_stdev_100[time],linestyle='-',lw=1,color = 'firebrick',label='r=100')
ax3.plot(x,rain_stdev_200[time],linestyle='-',lw=1,color = 'forestgreen',label='r=200 i='+str(time+5)+'hours a='+str("{:.2f}".format(10**rain_stdev_b[time]))) 
#ax3.plot(x,rain_stdev_500[time],linestyle='-',lw=1,color = 'mediumturquoise',label='r=500')  
    
ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.tick_params(axis='x', labelsize=19)
ax3.tick_params(axis='y', labelsize=19)
ax3.set_ylabel('$Width$ $of$ $95\%$ $CI$',fontsize=19);
ax3.set_xlabel('$Ensemble$ $size$ $(n)$',fontsize=19);
ax3.set_title('$Rain$ $stdev$',fontsize=19)

#################################################################################################
plt.tight_layout()
plt.savefig('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_plot_neighbourhood_eescomparison_stdev_extendedrun1_weakforcing_time_'+str(time)+'.pdf')

# plot q95 #

fig, (ax1, ax2, ax3)= plt.subplots(1, 3,figsize=[14,5])
plt.subplots_adjust(wspace = 0.5,hspace = 0)
plt.grid(True)
#################################################################################################

# do for only one time

ax1.plot(x,sf_wind_q95[time],linestyle='-',lw=4,alpha=0.5,color = 'dodgerblue',label='thicker=strong forcing 15') # strong forcing 
ax1.plot(x,sf_wind_q95_2[time],linestyle='-',lw=4,alpha=0.5,color = 'firebrick')
ax1.plot(x,sf_wind_q95_10[time],linestyle='-',lw=4,alpha=0.5,color = 'darkorange')
ax1.plot(x,sf_wind_q95_20[time],linestyle='-',lw=4,alpha=0.5,color = 'goldenrod')
ax1.plot(x,sf_wind_q95_50[time],linestyle='-',lw=4,alpha=0.5,color = 'chartreuse')
ax1.plot(x,sf_wind_q95_100[time],linestyle='-',lw=4,alpha=0.5,color = 'forestgreen')
ax1.plot(x,sf_wind_q95_200[time],linestyle='-',lw=4,alpha=0.5,color = 'rosybrown') 
#ax1.plot(x,sf_wind_q95_500[time],linestyle='-',lw=4,alpha=0.5,color = 'mediumturquoise')  

f = func(x, wind_q95_b[time])
ax1.fill_between(x,f+0.05*f,f-0.05*f,color='silver',alpha=0.5)
ax1.plot(x,wind_q95[time],linestyle='-',lw=1,color = 'dodgerblue',label='single gp')
ax1.plot(x,wind_q95_2[time],linestyle='-',lw=1,color = 'firebrick',label='r=2')
ax1.plot(x,wind_q95_10[time],linestyle='-',lw=1,color = 'darkorange',label='r=10')
ax1.plot(x,wind_q95_20[time],linestyle='-',lw=1,color = 'goldenrod',label='r=20')
ax1.plot(x,wind_q95_50[time],linestyle='-',lw=1,color = 'chartreuse',label='r=50')
ax1.plot(x,wind_q95_100[time],linestyle='-',lw=1,color = 'forestgreen',label='r=100')
ax1.plot(x,wind_q95_200[time],linestyle='-',lw=1,color = 'rosybrown',label='r=200 i='+str(time+5)+'hours a='+str("{:.2f}".format(10**wind_q95_b[time]))) 
#ax1.plot(x,wind_q95_500[time],linestyle='-',lw=1,color = 'mediumturquoise',label='r=500')  
    
ax1.legend()
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.tick_params(axis='x', labelsize=19)
ax1.tick_params(axis='y', labelsize=19)
ax1.set_ylabel('$Width$ $of$ $95\%$ $CI$',fontsize=19);
ax1.set_xlabel('$Ensemble$ $size$ $(n)$',fontsize=19);
ax1.set_title('$Wind$ $q95$',fontsize=19)

#################################################################################################

ax2.plot(x,sf_height_q95[time],linestyle='-',lw=4,alpha=0.5,color = 'dodgerblue',label='thicker=strong forcing 15') # strong forcing 
ax2.plot(x,sf_height_q95_2[time],linestyle='-',lw=4,alpha=0.5,color = 'firebrick')
ax2.plot(x,sf_height_q95_10[time],linestyle='-',lw=4,alpha=0.5,color = 'darkorange')
ax2.plot(x,sf_height_q95_20[time],linestyle='-',lw=4,alpha=0.5,color = 'goldenrod')
ax2.plot(x,sf_height_q95_50[time],linestyle='-',lw=4,alpha=0.5,color = 'chartreuse')
ax2.plot(x,sf_height_q95_100[time],linestyle='-',lw=4,alpha=0.5,color = 'forestgreen')
ax2.plot(x,sf_height_q95_200[time],linestyle='-',lw=4,alpha=0.5,color = 'rosybrown') 
#ax2.plot(x,sf_height_q95_500[time],linestyle='-',lw=4,alpha=0.5,color = 'mediumturquoise') 

f = func(x, height_q95_b[time])
ax2.fill_between(x,f+0.05*f,f-0.05*f,color='silver',alpha=0.5)
ax2.plot(x,height_q95[time],linestyle='-',lw=1,color = 'dodgerblue',label='single gp')
ax2.plot(x,height_q95_2[time],linestyle='-',lw=1,color = 'firebrick',label='r=2')
ax2.plot(x,height_q95_10[time],linestyle='-',lw=1,color = 'darkorange',label='r=10')
ax2.plot(x,height_q95_20[time],linestyle='-',lw=1,color = 'goldenrod',label='r=20')
ax2.plot(x,height_q95_50[time],linestyle='-',lw=1,color = 'chartreuse',label='r=50')
ax2.plot(x,height_q95_100[time],linestyle='-',lw=1,color = 'forestgreen',label='r=100')
ax2.plot(x,height_q95_200[time],linestyle='-',lw=1,color = 'rosybrown',label='r=200 i='+str(time+5)+'hours a='+str("{:.2f}".format(10**height_q95_b[time]))) 
#ax2.plot(x,height_q95_500[time],linestyle='-',lw=1,color = 'mediumturquoise',label='r=500')  
    
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.tick_params(axis='x', labelsize=19)
ax2.tick_params(axis='y', labelsize=19)
ax2.set_ylabel('$Width$ $of$ $95\%$ $CI$',fontsize=19);
ax2.set_xlabel('$Ensemble$ $size$ $(n)$',fontsize=19);
ax2.set_title('$Height$ $q95$',fontsize=19)

#################################################################################################

ax3.plot(x,sf_rain_q95[time],linestyle='-',lw=4,alpha=0.5,color = 'dodgerblue',label='thicker=strong forcing 15') # strong forcing 
ax3.plot(x,sf_rain_q95_2[time],linestyle='-',lw=4,alpha=0.5,color = 'firebrick')
ax3.plot(x,sf_rain_q95_10[time],linestyle='-',lw=4,alpha=0.5,color = 'darkorange')
ax3.plot(x,sf_rain_q95_20[time],linestyle='-',lw=4,alpha=0.5,color = 'goldenrod')
ax3.plot(x,sf_rain_q95_50[time],linestyle='-',lw=4,alpha=0.5,color = 'chartreuse')
ax3.plot(x,sf_rain_q95_100[time],linestyle='-',lw=4,alpha=0.5,color = 'forestgreen')
ax3.plot(x,sf_rain_q95_200[time],linestyle='-',lw=4,alpha=0.5,color = 'rosybrown') 
#ax3.plot(x,sf_rain_q95_500[time],linestyle='-',lw=4,alpha=0.5,color = 'mediumturquoise') 

f = func(x, rain_q95_b[time])
ax3.fill_between(x,f+0.05*f,f-0.05*f,color='silver',alpha=0.5)
ax3.plot(x,rain_q95[time],linestyle='-',lw=1,color = 'dodgerblue',label='single gp')
ax3.plot(x,rain_q95_2[time],linestyle='-',lw=1,color = 'firebrick',label='r=2')
ax3.plot(x,rain_q95_10[time],linestyle='-',lw=1,color = 'darkorange',label='r=10')
ax3.plot(x,rain_q95_20[time],linestyle='-',lw=1,color = 'goldenrod',label='r=20')
ax3.plot(x,rain_q95_50[time],linestyle='-',lw=1,color = 'chartreuse',label='r=50')
ax3.plot(x,rain_q95_100[time],linestyle='-',lw=1,color = 'forestgreen',label='r=100')
ax3.plot(x,rain_q95_200[time],linestyle='-',lw=1,color = 'rosybrown',label='r=200 i='+str(time+5)+'hours a='+str("{:.2f}".format(10**rain_q95_b[time]))) 
#ax3.plot(x,rain_q95_500[time],linestyle='-',lw=1,color = 'mediumturquoise',label='r=500')  
    
ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.tick_params(axis='x', labelsize=19)
ax3.tick_params(axis='y', labelsize=19)
ax3.set_ylabel('$Width$ $of$ $95\%$ $CI$',fontsize=19);
ax3.set_xlabel('$Ensemble$ $size$ $(n)$',fontsize=19);
ax3.set_title('$Rain$ $q95$',fontsize=19)

#################################################################################################
plt.tight_layout()
plt.savefig('/scratch/k/K.Tempest/Model_extension/weak_forcing_1/convergence_data/convergence_plot_neighbourhood_eescomparison_q95_extendedrun1_weakforcing_time_'+str(time)+'.pdf')

