import numpy as np
import os
import argparse
import numpy.ma as ma
from netCDF4 import Dataset
import pickle

################################################

# creates distribution for a specific grid point and for every timestep. 
# distribution is of full ensemble. 
# save file with dimensions: []

################################################

parser = argparse.ArgumentParser(description='Sorts distribution data. Input variables please.')
parser.add_argument('--grid_point',required=True, type=int, help='grid point')
args = parser.parse_args()
grid_point = args.grid_point
    
#########################################################
# variables: 


#########################################################

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

# phi_c data (only compute once as it is the same for every gridpoint across the domain)

data_phi_c_0 = np.concatenate((ma.getdata(SWM00['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM600['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM1200['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM1800['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM2400['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM3000['/free_run/truth_phi_c'][grid_point,:,:,:])),axis=2)
data_phi_c_1 = np.concatenate((ma.getdata(SWM01['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM601['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM1201['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM1801['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM2401['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM3001['/free_run/truth_phi_c'][grid_point,:,:,:])),axis=2)
data_phi_c_2 = np.concatenate((ma.getdata(SWM02['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM602['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM1202['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM1802['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM2402['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM3002['/free_run/truth_phi_c'][grid_point,:,:,:])),axis=2)
data_phi_c_3 = np.concatenate((ma.getdata(SWM03['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM603['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM1203['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM1803['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM2403['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM3003['/free_run/truth_phi_c'][grid_point,:,:,:])),axis=2)
data_phi_c_4 = np.concatenate((ma.getdata(SWM04['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM604['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM1204['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM1804['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM2404['/free_run/truth_phi_c'][grid_point,:,:,:]),ma.getdata(SWM3004['/free_run/truth_phi_c'][grid_point,:,:,:])),axis=2)

data_phi_c = np.concatenate((data_phi_c_0,data_phi_c_1,data_phi_c_2,data_phi_c_3,data_phi_c_4),axis=1)

data_phi_c = data_phi_c[1,:,:]

# phi data 

data_phi_0 = np.concatenate((ma.getdata(SWM00['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM600['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM1200['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM1800['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM2400['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM3000['/free_run/truth_phi'][grid_point,:,:])),axis=1)
data_phi_1 = np.concatenate((ma.getdata(SWM01['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM601['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM1201['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM1801['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM2401['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM3001['/free_run/truth_phi'][grid_point,:,:])),axis=1)
data_phi_2 = np.concatenate((ma.getdata(SWM02['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM602['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM1202['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM1802['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM2402['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM3002['/free_run/truth_phi'][grid_point,:,:])),axis=1)
data_phi_3 = np.concatenate((ma.getdata(SWM03['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM603['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM1203['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM1803['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM2403['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM3003['/free_run/truth_phi'][grid_point,:,:])),axis=1)
data_phi_4 = np.concatenate((ma.getdata(SWM04['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM604['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM1204['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM1804['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM2404['/free_run/truth_phi'][grid_point,:,:]),ma.getdata(SWM3004['/free_run/truth_phi'][grid_point,:,:])),axis=1)

data_phi = np.concatenate((data_phi_0,data_phi_1,data_phi_2,data_phi_3,data_phi_4),axis=0)

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
data_w = 0.5*(np.sum([a[1:n+1,:,:]+a[0:n,:,:]],axis=0))
data_w = data_w[grid_point,:,:]

# free up space 

data_w_0 = data_w_1 = data_w_2 = data_w_3 = data_w_4 = 0

# height data

data_h_0 = np.concatenate((ma.getdata(SWM00['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM600['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM1200['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM1800['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM2400['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM3000['/free_run/truth'][grid_point + 1000,:,:])),axis=1)
data_h_1 = np.concatenate((ma.getdata(SWM01['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM601['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM1201['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM1801['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM2401['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM3001['/free_run/truth'][grid_point + 1000,:,:])),axis=1)
data_h_2 = np.concatenate((ma.getdata(SWM02['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM602['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM1202['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM1802['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM2402['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM3002['/free_run/truth'][grid_point + 1000,:,:])),axis=1)
data_h_3 = np.concatenate((ma.getdata(SWM03['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM603['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM1203['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM1803['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM2403['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM3003['/free_run/truth'][grid_point + 1000,:,:])),axis=1)
data_h_4 = np.concatenate((ma.getdata(SWM04['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM604['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM1204['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM1804['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM2404['/free_run/truth'][grid_point + 1000,:,:]),ma.getdata(SWM3004['/free_run/truth'][grid_point + 1000,:,:])),axis=1)

data_h = np.concatenate((data_h_0,data_h_1,data_h_2,data_h_3,data_h_4),axis=0)

# rain data

data_r_0 = np.concatenate((ma.getdata(SWM00['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM600['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM1200['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM1800['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM2400['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM3000['/free_run/truth'][grid_point + 2000,:,:])),axis=1)
data_r_1 = np.concatenate((ma.getdata(SWM01['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM601['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM1201['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM1801['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM2401['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM3001['/free_run/truth'][grid_point + 2000,:,:])),axis=1)
data_r_2 = np.concatenate((ma.getdata(SWM02['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM602['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM1202['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM1802['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM2402['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM3002['/free_run/truth'][grid_point + 2000,:,:])),axis=1)
data_r_3 = np.concatenate((ma.getdata(SWM03['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM603['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM1203['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM1803['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM2403['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM3003['/free_run/truth'][grid_point + 2000,:,:])),axis=1)
data_r_4 = np.concatenate((ma.getdata(SWM04['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM604['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM1204['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM1804['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM2404['/free_run/truth'][grid_point + 2000,:,:]),ma.getdata(SWM3004['/free_run/truth'][grid_point + 2000,:,:])),axis=1)

data_r = np.concatenate((data_r_0,data_r_1,data_r_2,data_r_3,data_r_4),axis=0)

# save data
with open('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/distribution_data/distribution_data_weakforcing275000_gp_'+str(grid_point)+'_phi_c','wb') as f:
    pickle.dump(data_phi_c,f)

with open('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/distribution_data/distribution_data_weakforcing275000_gp_'+str(grid_point)+'_phi','wb') as f:
    pickle.dump(data_phi,f)

with open('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/distribution_data/distribution_data_weakforcing275000_gp_'+str(grid_point)+'_wind','wb') as f:
    pickle.dump(data_w,f)

with open('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/distribution_data/distribution_data_weakforcing275000_gp_'+str(grid_point)+'_height','wb') as f:
    pickle.dump(data_h,f)

with open('/project/meteo/scratch/K.Tempest/Model_extension_large/weak_forcing_27_5000/distribution_data/distribution_data_weakforcing275000_gp_'+str(grid_point)+'_rain','wb') as f:
    pickle.dump(data_r,f)


