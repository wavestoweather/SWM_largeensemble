from msw_model_expanded import *
from DA_2019 import *
from constants import *
from Plotting import *
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import os
import argparse

################################################
# get start_time
parser = argparse.ArgumentParser(description='Running DA and/or free run of SWM. Input start time please.')
parser.add_argument('--start_time',required=False, type=int, help='start time')
parser.add_argument('--batch_number',required=False, type=int, help='batch_number')
parser.add_argument('--members_per_batch',required=False, type=int, help='members_per_batch')
parser.add_argument('--DA',required=True, type=int, help='If want DA, use 1')
args = parser.parse_args()
start_time = args.start_time
batch_number = args.batch_number
members_per_batch = args.members_per_batch
DA = args.DA

#################################################
# DA

if DA==1:
    DA_start = datetime.now()
    print('starting DA') # DA initialisation
    DA_init_vis_dt__(EnKF)
    print('finished DA in time of = '+str(DA_start-datetime.now()))

model_start_time = datetime.now()
# free model run
model_initss__(start_time=start_time,batch_number=batch_number,members_per_batch=members_per_batch)
print('done model run with batch number = '+str(batch_number)+ 'and members per batch of ' +str(members_per_batch)+' in time '+str(model_start_time-datetime.now()))
