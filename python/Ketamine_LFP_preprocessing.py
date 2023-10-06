# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from utilities_ketamine_analysis_v8 import *
from utils_signal_processing import *
from utils_plotting import *
from utils_general import *


#%%

sess = 2 # session number 

binFullPath = r'C:\Users\fentonlab\Desktop\Gino\LFPs'
HPC_path_file = os.path.join(r'C:\Users\fentonlab\Desktop\Gino\LFPs','HPC_lfp_paths.file')
PFC_path_file = os.path.join(r'C:\Users\fentonlab\Desktop\Gino\LFPs','PFC_lfp_paths.file')

#%%

# =============================================================================
# Load LFP and speed - upsample speed 
# =============================================================================
# Load Lfp and speed data for a specific recording and brain area 
Lfp, speed = load_data(binFullPath,HPC_path_file,PFC_path_file,"HPC",sess)

# Upsample Speed recording
speed_up = upsample_speed(speed, Lfp, sess, LFP_rate = 2500, speed_rate = 100)


# =============================================================================
# PLOTTING -- sanity checks / data visualization 
# =============================================================================

# plot speed
plot_speed(speed_up,sess,1,100,10,fs_lfp=2500,fs_beh=100)
# plot speed histogram 
plot_speed_histo_logscale(speed,sess)
plot_speed_histo_regimes(speed_up, th_low = 50,th_mid = 100)

# plot Lfp
plot_lfp_two_channels(Lfp,10,36,0,100,10,N=2500)
plot_lfp_various_channels(Lfp,1,9,10,100,3,3,10,N=2500)
 
 
#%%

# 2. Split speed and Lfp into injection periods (epochs): baseline, low, mid, and high injection

Lfp_B, Lfp_L, Lfp_M, Lfp_H, speed_B, speed_L, speed_M, speed_H = split_into_epochs(Lfp,speed_up,N=2500)
 
#%%
#### Create list to store Lfp for each epoch: (channel, minute, n trial, trial data )
 
nch = lfp.shape[1]
lfp_B_epoch_low_s = [[] for ch in range(nch)]
lfp_L_epoch_low_s = [[] for ch in range(nch)]
lfp_M_epoch_low_s = [[] for ch in range(nch)]
lfp_H_epoch_low_s = [[] for ch in range(nch)]


lfp_B_epoch_high_s = [[] for ch in range(nch)]
lfp_L_epoch_high_s = [[] for ch in range(nch)]
lfp_M_epoch_high_s = [[] for ch in range(nch)]
lfp_H_epoch_high_s = [[] for ch in range(nch)]





#### 4. Filter Lfp in each epoch

filter_lfp_in_each_epoch(Lfp_B_min,Lfp_L_min,Lfp_M_min,Lfp_H_min,gain)






 
 
 