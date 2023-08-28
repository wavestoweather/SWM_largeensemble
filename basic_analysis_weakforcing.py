import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import stats
from datetime import datetime
from numpy import ndarray
from netCDF4 import Dataset
import numpy.ma as ma
import pickle
from datetime import datetime
import seaborn as sns
import matplotlib.animation as anim
import scipy.stats as ss
import pandas as pd
from scipy.stats import kurtosis, skew, ks_2samp, kstest
import random
from sklearn.utils import shuffle

sns.set(color_codes = True)
sns.set_style('whitegrid')

# import data 

SWM00 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM0_0_1000.nc','r') 
SWM01 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM0_1_1000.nc','r') 
SWM02 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM0_2_1000.nc','r') 
SWM03 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM0_3_1000.nc','r') 
SWM04 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM0_4_1000.nc','r') 

SWM600 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM60_0_1000.nc','r') 
SWM601 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM60_1_1000.nc','r') 
SWM602 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM60_2_1000.nc','r') 
SWM603 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM60_3_1000.nc','r') 
SWM604 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM60_4_1000.nc','r') 

SWM1200 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM120_0_1000.nc','r') 
SWM1201 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM120_1_1000.nc','r') 
SWM1202 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM120_2_1000.nc','r') 
SWM1203 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM120_3_1000.nc','r') 
SWM1204 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM120_4_1000.nc','r') 

SWM1800 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM180_0_1000.nc','r') 
SWM1801 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM180_1_1000.nc','r') 
SWM1802 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM180_2_1000.nc','r') 
SWM1803 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM180_3_1000.nc','r') 
SWM1804 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM180_4_1000.nc','r') 

SWM2400 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM240_0_1000.nc','r') 
SWM2401 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM240_1_1000.nc','r') 
SWM2402 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM240_2_1000.nc','r') 
SWM2403 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM240_3_1000.nc','r') 
SWM2404 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM240_4_1000.nc','r') 

SWM3000 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM300_0_1000.nc','r') 
SWM3001 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM300_1_1000.nc','r') 
SWM3002 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM300_2_1000.nc','r') 
SWM3003 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM300_3_1000.nc','r') 
SWM3004 = Dataset('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/SWM300_4_1000.nc','r') 

# for individual groups first
data_phi_c_0 = np.concatenate((ma.getdata(SWM00['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM600['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM1200['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM1800['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM2400['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM3000['/free_run/truth_phi_c'][:,:,:,:])),axis=3)
data_phi_0 = np.concatenate((ma.getdata(SWM00['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM600['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM1200['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM1800['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM2400['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM3000['/free_run/truth_phi'][:,:,:])),axis=2) 

data_phi_c_1 = np.concatenate((ma.getdata(SWM01['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM601['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM1201['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM1801['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM2401['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM3001['/free_run/truth_phi_c'][:,:,:,:])),axis=3)
data_phi_1 = np.concatenate((ma.getdata(SWM01['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM601['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM1201['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM1801['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM2401['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM3001['/free_run/truth_phi'][:,:,:])),axis=2)                        

data_phi_c_2 = np.concatenate((ma.getdata(SWM02['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM602['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM1202['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM1802['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM2402['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM3002['/free_run/truth_phi_c'][:,:,:,:])),axis=3)
data_phi_2 = np.concatenate((ma.getdata(SWM02['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM602['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM1202['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM1802['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM2402['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM3002['/free_run/truth_phi'][:,:,:])),axis=2)                        

data_phi_c_3 = np.concatenate((ma.getdata(SWM03['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM603['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM1203['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM1803['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM2403['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM3003['/free_run/truth_phi_c'][:,:,:,:])),axis=3)
data_phi_3 = np.concatenate((ma.getdata(SWM03['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM603['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM1203['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM1803['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM2403['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM3003['/free_run/truth_phi'][:,:,:])),axis=2)                        

data_phi_c_4 = np.concatenate((ma.getdata(SWM04['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM604['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM1204['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM1804['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM2404['/free_run/truth_phi_c'][:,:,:,:]),ma.getdata(SWM3004['/free_run/truth_phi_c'][:,:,:,:])),axis=3)
data_phi_4 = np.concatenate((ma.getdata(SWM04['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM604['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM1204['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM1804['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM2404['/free_run/truth_phi'][:,:,:]),ma.getdata(SWM3004['/free_run/truth_phi'][:,:,:])),axis=2)                        

data_phi_c = np.concatenate((data_phi_c_0,data_phi_c_1,data_phi_c_2,data_phi_c_3,data_phi_c_4),axis=2)
data_phi = np.concatenate((data_phi_0,data_phi_1,data_phi_2,data_phi_3,data_phi_4),axis=1)

data_phi_c_0 = data_phi_c_1 = data_phi_c_2 = data_phi_c_3 = data_phi_c_4 = 0
data_phi_0 = data_phi_1 = data_phi_2 = data_phi_3 = data_phi_4 = 0

data_phi_c = data_phi_c[:,1,:,:]

data_av_phi = np.mean(data_phi,axis=0)
data_av_phi_c = np.mean(data_phi_c[:,:,:],axis=0)

x = np.linspace(0,23,360) # when 24 hours there are 360 datapoints 
x_24 = np.arange(1,25,1)
n = 1000

print('loaded phi data')

# distribution data
# wind 
data_w_0 = np.concatenate((ma.getdata(SWM00['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM600['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1200['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1800['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM2400['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM3000['/free_run/truth'][0:1000,:,:])),axis=2)
data_w_1 = np.concatenate((ma.getdata(SWM01['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM601['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1201['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1801['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM2401['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM3001['/free_run/truth'][0:1000,:,:])),axis=2)
data_w_2 = np.concatenate((ma.getdata(SWM02['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM602['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1202['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1802['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM2402['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM3002['/free_run/truth'][0:1000,:,:])),axis=2)
data_w_3 = np.concatenate((ma.getdata(SWM03['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM603['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1203['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1803['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM2403['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM3003['/free_run/truth'][0:1000,:,:])),axis=2)
data_w_4 = np.concatenate((ma.getdata(SWM04['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM604['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1204['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM1804['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM2404['/free_run/truth'][0:1000,:,:]),ma.getdata(SWM3004['/free_run/truth'][0:1000,:,:])),axis=2)

data_w = np.concatenate((data_w_0,data_w_1,data_w_2,data_w_3,data_w_4),axis=1)

# interpolate wind data 
n = 1000
g = 10
a = np.insert(data_w, n, data_w[0,:,:],axis=0)
data_w = 0.5*(np.sum([a[1:n+1,:,:]+a[0:n,:,:]],axis=0))

data_w_mean = np.mean(data_w[:,:,:],axis=0)

# to free up space 
data_w_0 = data_w_1 = data_w_2 = data_w_3 = data_w_4 = 0

# height
dp = np.arange(1000,2000)

data_h_0 = np.concatenate((ma.getdata(SWM00['/free_run/truth'][dp,:,:]),ma.getdata(SWM600['/free_run/truth'][dp,:,:]),ma.getdata(SWM1200['/free_run/truth'][dp,:,:]),ma.getdata(SWM1800['/free_run/truth'][dp,:,:]),ma.getdata(SWM2400['/free_run/truth'][dp,:,:]),ma.getdata(SWM3000['/free_run/truth'][dp,:,:])),axis=2)
data_h_1 = np.concatenate((ma.getdata(SWM01['/free_run/truth'][dp,:,:]),ma.getdata(SWM601['/free_run/truth'][dp,:,:]),ma.getdata(SWM1201['/free_run/truth'][dp,:,:]),ma.getdata(SWM1801['/free_run/truth'][dp,:,:]),ma.getdata(SWM2401['/free_run/truth'][dp,:,:]),ma.getdata(SWM3001['/free_run/truth'][dp,:,:])),axis=2)
data_h_2 = np.concatenate((ma.getdata(SWM02['/free_run/truth'][dp,:,:]),ma.getdata(SWM602['/free_run/truth'][dp,:,:]),ma.getdata(SWM1202['/free_run/truth'][dp,:,:]),ma.getdata(SWM1802['/free_run/truth'][dp,:,:]),ma.getdata(SWM2402['/free_run/truth'][dp,:,:]),ma.getdata(SWM3002['/free_run/truth'][dp,:,:])),axis=2)
data_h_3 = np.concatenate((ma.getdata(SWM03['/free_run/truth'][dp,:,:]),ma.getdata(SWM603['/free_run/truth'][dp,:,:]),ma.getdata(SWM1203['/free_run/truth'][dp,:,:]),ma.getdata(SWM1803['/free_run/truth'][dp,:,:]),ma.getdata(SWM2403['/free_run/truth'][dp,:,:]),ma.getdata(SWM3003['/free_run/truth'][dp,:,:])),axis=2)
data_h_4 = np.concatenate((ma.getdata(SWM04['/free_run/truth'][dp,:,:]),ma.getdata(SWM604['/free_run/truth'][dp,:,:]),ma.getdata(SWM1204['/free_run/truth'][dp,:,:]),ma.getdata(SWM1804['/free_run/truth'][dp,:,:]),ma.getdata(SWM2404['/free_run/truth'][dp,:,:]),ma.getdata(SWM3004['/free_run/truth'][dp,:,:])),axis=2)

data_h = np.concatenate((data_h_0,data_h_1,data_h_2,data_h_3,data_h_4),axis=1)

data_h_mean = np.mean(data_h[:,:,:],axis=0)

data_h_mean_hm = np.mean(data_h[:,:,:],axis=1)
print('loaded height data')
# rain 
dp = np.arange(2000,3000)

data_r_0 = np.concatenate((ma.getdata(SWM00['/free_run/truth'][dp,:,:]),ma.getdata(SWM600['/free_run/truth'][dp,:,:]),ma.getdata(SWM1200['/free_run/truth'][dp,:,:]),ma.getdata(SWM1800['/free_run/truth'][dp,:,:]),ma.getdata(SWM2400['/free_run/truth'][dp,:,:]),ma.getdata(SWM3000['/free_run/truth'][dp,:,:])),axis=2)
data_r_1 = np.concatenate((ma.getdata(SWM01['/free_run/truth'][dp,:,:]),ma.getdata(SWM601['/free_run/truth'][dp,:,:]),ma.getdata(SWM1201['/free_run/truth'][dp,:,:]),ma.getdata(SWM1801['/free_run/truth'][dp,:,:]),ma.getdata(SWM2401['/free_run/truth'][dp,:,:]),ma.getdata(SWM3001['/free_run/truth'][dp,:,:])),axis=2)
data_r_2 = np.concatenate((ma.getdata(SWM02['/free_run/truth'][dp,:,:]),ma.getdata(SWM602['/free_run/truth'][dp,:,:]),ma.getdata(SWM1202['/free_run/truth'][dp,:,:]),ma.getdata(SWM1802['/free_run/truth'][dp,:,:]),ma.getdata(SWM2402['/free_run/truth'][dp,:,:]),ma.getdata(SWM3002['/free_run/truth'][dp,:,:])),axis=2)
data_r_3 = np.concatenate((ma.getdata(SWM03['/free_run/truth'][dp,:,:]),ma.getdata(SWM603['/free_run/truth'][dp,:,:]),ma.getdata(SWM1203['/free_run/truth'][dp,:,:]),ma.getdata(SWM1803['/free_run/truth'][dp,:,:]),ma.getdata(SWM2403['/free_run/truth'][dp,:,:]),ma.getdata(SWM3003['/free_run/truth'][dp,:,:])),axis=2)
data_r_4 = np.concatenate((ma.getdata(SWM04['/free_run/truth'][dp,:,:]),ma.getdata(SWM604['/free_run/truth'][dp,:,:]),ma.getdata(SWM1204['/free_run/truth'][dp,:,:]),ma.getdata(SWM1804['/free_run/truth'][dp,:,:]),ma.getdata(SWM2404['/free_run/truth'][dp,:,:]),ma.getdata(SWM3004['/free_run/truth'][dp,:,:])),axis=2)

data_r = np.concatenate((data_r_0,data_r_1,data_r_2,data_r_3,data_r_4),axis=1)

data_r_0 = data_r_1 = data_r_2 = data_r_3 = data_r_4 = 0 # free up space

data_r_mean = np.mean(data_r[:,:,:],axis=0) #average over grid points

data_r_mean_low = np.quantile(data_r_mean,q=0.025,axis=0)
data_r_mean_high = np.quantile(data_r_mean,q=0.975,axis=0)

print('loaded all data now')

# get CAPE confidence interval:

CAPE_low = np.quantile(g*data_h_mean[:,:]-data_av_phi_c[:,:],q=0.025,axis=0)
CAPE_high = np.quantile(g*data_h_mean[:,:]-data_av_phi_c[:,:],q=0.975,axis=0)

# plot rain and CAPE curve # for an averaged domain (new definition)
g = 10

fig = plt.figure(figsize=[14,6])

# precipitation

# for i in range(100):
#     plt.plot(x,data_r_mean[(i+1)*50-1,:]*2000,'.',color='darkblue',alpha=0.1)
plt.fill_between(x,data_r_mean_high*2000,data_r_mean_low*2000,color='dodgerblue',alpha=0.3)
plt.plot(x,np.mean(data_r_mean[:,:],axis=0)*2000,'.',color='dodgerblue',label='rain')

# CAPE

# for i in range(100):
#     plt.plot(x,g*data_h_mean[(i+1)*50-1,:]-data_av_phi_c[(i+1)*50-1,:],'darkred',alpha=0.1)
plt.fill_between(x,CAPE_high,CAPE_low,color='darkred',alpha=0.3)
plt.plot(x,np.mean(g*data_h_mean[:,:]-data_av_phi_c[:,:],axis=0),'darkred',label='CAPE')

plt.xlabel('$Hour$',fontsize=16)
plt.ylabel('')
plt.title('$CAPE$ $and$ $precipitation$ $(x2000)$',fontsize=16)
plt.legend(fontsize=16);
plt.tight_layout()

plt.savefig('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/general/CAPE_and_rain_NEWCAPE.pdf')

# HM Diagram 
# single
plt.cla()

fig, ax = plt.subplots(figsize=(14, 8))
#levels = [37.9,38,38.1,38.2,38.3,38.4,38.5,38.6,38.7,38.8,38.9]

total = data_h[:,24,:] # 1 ensemble member
CS = plt.contourf(total,cmap=plt.cm.YlGn)#,levels=levels);

ci = fig.colorbar(CS,orientation="horizontal", pad=0.14)
ci.ax.tick_params(labelsize=18)

plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=18)

ax.set_xticks(np.linspace(14,359,24));
ax.set_xticklabels(np.arange(1,25));

plt.xlabel('$time$ $(hour)$',fontsize=19);
plt.ylabel('$domain$ $(km)$',fontsize=19);
ci.ax.set_xlabel('$height$ $(m)$', rotation=0,fontsize=19)

plt.tight_layout()

plt.savefig('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/general/Height_single_HM.pdf')

# averaged HM

# mean over ensemble member - height #
plt.cla()
fig, ax = plt.subplots(figsize=(14, 8))

total = data_h_mean_hm[:,:] # 1 ensemble member
CS = plt.contourf(total,cmap=plt.cm.YlGn)#,levels=levels);

ci = fig.colorbar(CS,orientation="horizontal", pad=0.14)
ci.ax.tick_params(labelsize=18)

plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=18)

ax.set_xticks(np.linspace(14,359,24));
ax.set_xticklabels(np.arange(1,25));

plt.xlabel('$time$ $(hour)$',fontsize=19);
plt.ylabel('$domain$ $(km)$',fontsize=19);
ci.ax.set_xlabel('$height$ $(m)$', rotation=0,fontsize=19)

plt.tight_layout()

plt.savefig('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/general/Height_averaged_HM.pdf')

# calculate convective timescale
# constants
gamma_2 = -2000
beta = 0.1
n = 1000
dx = 500
g = 10

# calculate conv_av
data_w_extra = np.insert(data_w, 0, data_w[n-1,:,:],axis=0)
data_w_extra_2 = np.insert(data_w_extra, n+1, data_w_extra[1,:,:],axis=0)
array = np.where(data_h<38.4,0,(data_w_extra_2[2:n+2,:,:] - data_w_extra_2[1:n+1,:,:])/dx)
array = np.where(array[:,:]>0,0,array[:,:])
con_av = np.mean(array[:,:],axis=0)
con_av = np.nan_to_num(con_av)

diffCAPE = -gamma_2 * beta * array[:,:,:]*3600
diffCAPEdomav = np.mean(diffCAPE,axis=0)
CAPE = np.mean(g*data_h_mean[:,:]-data_av_phi_c[:,:],axis=0)

# Way 1: average over ensemble before calculating fraction
plt.cla()
fig = plt.figure(figsize=[10,4])
plt.plot(x,-CAPE/np.mean(diffCAPEdomav,axis=0),'.',color='g',label='T_c');
plt.title('tau_c equation');
plt.legend();
plt.xlabel('hours')
plt.ylabel('T_c')
plt.ylim(0,40)

plt.savefig('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/general/Convective_timescale.pdf')

plt.cla()
fig, ax1 = plt.subplots(figsize=[9,7])

color = 'tab:red'

ax1.set_xlabel('$free$ $run$ $(hours)$',fontsize=20)
ax1.set_ylabel('$CAPE$ $(m^2s^{-2})$ $and$ $rain$ $(x2000)$',fontsize=20)

ax1.fill_between(x,data_r_mean_high*2000,data_r_mean_low*2000,color='dodgerblue',alpha=0.3)
ax1.plot(x,np.mean(data_r_mean[:,:],axis=0)*2000,color='dodgerblue',label='rain')

# CAPE
ax1.fill_between(x,CAPE_high,CAPE_low,color='darkred',alpha=0.3)
ax1.plot(x,np.mean(g*data_h_mean[:,:]-data_av_phi_c[:,:],axis=0),'darkred',label='CAPE')
ax1.set_ylim(0,0.73)

CAPE = -CAPE/np.mean(diffCAPEdomav,axis=0) # prepare simple moving average
CAPE_av = np.empty(24)
for i in range(24):
    CAPE_av[i] = np.mean(CAPE[(i*15):((i+1)*15)])

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:green'
ax2.set_ylabel('$Convective$ $timescale$ $(hours)$', color=color,fontsize=20)  # we already handled the x-label with ax1
ax2.plot(x_24, CAPE_av, color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(0,100)

ax1.legend(loc=1,fontsize=14)

plt.savefig('/scratch/k/K.Tempest/Model_extension_large/weak_forcing_27_5000/general/CAPE_prec_convective_timescale_paper.pdf')

