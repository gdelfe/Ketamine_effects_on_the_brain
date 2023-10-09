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


# =============================================================================
# Load LFP and speed - upsample speed 
# =============================================================================
# Load Lfp and speed data for a specific recording and brain area 
Lfp, speed, gain, rec = load_data(binFullPath,HPC_path_file,PFC_path_file,"HPC",sess)

# Upsample Speed recording
speed_up = upsample_speed(speed, Lfp, sess, LFP_rate = 2500, speed_rate = 100)

# =============================================================================
# PLOTTING -- sanity checks / data visualization 
# =============================================================================

# # plot speed
# plot_speed(speed_up,sess,1,100,10,fs_lfp=2500,fs_beh=100)
# # plot speed histogram 
# plot_speed_histo_logscale(speed,sess)
# plot_speed_histo_regimes(speed_up, th_low = 50,th_mid = 100)

# # plot Lfp
# plot_lfp_two_channels(Lfp,10,36,0,100,10,N=2500)
# plot_lfp_various_channels(Lfp,1,9,10,100,3,3,10,N=2500)
 
 
#%%

# Split speed and Lfp into injection periods (epochs): baseline, low, mid, and high injection

Lfp_B, Lfp_L, Lfp_M, Lfp_H, speed_B, speed_L, speed_M, speed_H = split_into_epochs(Lfp,speed_up,N=2500)
 

# Create list to store Lfp for each epoch: (channel, minute, n trial, trial data )
 
nch = Lfp_B.shape[1]
# low speed
lfp_B_ep_low_s = [[] for ch in range(nch)]
lfp_L_ep_low_s = [[] for ch in range(nch)]
lfp_M_ep_low_s = [[] for ch in range(nch)]
lfp_H_ep_low_s = [[] for ch in range(nch)]
# high speed
lfp_B_ep_high_s = [[] for ch in range(nch)]
lfp_L_ep_high_s = [[] for ch in range(nch)]
lfp_M_ep_high_s = [[] for ch in range(nch)]
lfp_H_ep_high_s = [[] for ch in range(nch)]

# =============================================================================
# SELECT ONE MINUTE DATA 
# =============================================================================

# Select 1 min data 
current_min = 0
Lfp_B_min, Lfp_L_min, Lfp_M_min, Lfp_H_min, speed_B_min, speed_L_min, speed_M_min, speed_H_min = select_1min_data(Lfp_B, Lfp_L, Lfp_M, Lfp_H, speed_B, speed_L, speed_M, speed_H, current_min, N=2500)

# =============================================================================
# Filter LFP (band pass)
# =============================================================================

# Filter Lfp in each epoch
lfp_filt_B, lfp_filt_L, lfp_filt_M, lfp_filt_H = filter_lfp_in_each_epoch(Lfp_B_min,Lfp_L_min,Lfp_M_min,Lfp_H_min,gain)

# plot_filtered_lfp(lfp_filt_B,0,10,1,36, 2500)

# Decimate Lfp and speed (subsample)
lfp_dec_B, lfp_dec_L, lfp_dec_M, lfp_dec_H, speed_dec_B, speed_dec_L, speed_dec_M, speed_dec_H = decimate_lfp_and_speed(lfp_filt_B,lfp_filt_L,lfp_filt_M,lfp_filt_H,speed_B_min,speed_L_min,speed_M_min,speed_H_min)


# =============================================================================
# MASKING SPEED AND LFP ARTIFACTS 
# =============================================================================

tot_mask_B_low_s, tot_mask_L_low_s, tot_mask_M_low_s, tot_mask_H_low_s, tot_mask_B_high_s, tot_mask_L_high_s, tot_mask_M_high_s,tot_mask_H_high_s = make_speed_and_lfp_maks(lfp_dec_B,lfp_dec_L, lfp_dec_M, lfp_dec_H, speed_dec_B, speed_dec_L, speed_dec_M, speed_dec_H, win = 1250, th = 30)

# =============================================================================
#  Reshape Lfp into: trial number, trial length, channels
# =============================================================================

 # # reshape Lfp in trial x time window x channels
win = 1250
LfpRB = lfp_dec_B.reshape(int(lfp_dec_B.shape[0]/win),-1,lfp_dec_B.shape[1]) # reshape Lfp: trial x win length x channel
LfpRL = lfp_dec_L.reshape(int(lfp_dec_L.shape[0]/win),-1,lfp_dec_L.shape[1]) # reshape Lfp: trial x win length x channel
LfpRM = lfp_dec_M.reshape(int(lfp_dec_M.shape[0]/win),-1,lfp_dec_M.shape[1]) # reshape Lfp: trial x win length x channel
LfpRH = lfp_dec_H.reshape(int(lfp_dec_H.shape[0]/win),-1,lfp_dec_H.shape[1]) # reshape Lfp: trial x win length x channel

print('reshaped Lfp, ', LfpRB.shape, LfpRL.shape, LfpRM.shape, LfpRH.shape)

#%%

# =============================================================================
# Keep good trials only and stack them into a 4D array
# =============================================================================

# =============================================================================
# LOW SPEED 

# keep only good trial for low speed
lfp_B_low_s_list, lfp_L_low_s_list, lfp_M_low_s_list, lfp_H_low_s_list = keep_only_good_trials(LfpRB, LfpRL, LfpRM, LfpRH, tot_mask_B_low_s, tot_mask_L_low_s, tot_mask_M_low_s, tot_mask_H_low_s, "low speed")
# stack 1 min Lfp for low speed trials into a 4D list/array: nch, min id, id trial, length trial,
lfp_B_ep_low_s, lfp_L_ep_low_s, lfp_M_ep_low_s, lfp_H_ep_low_s = stack_lfp_1min(lfp_B_ep_low_s, lfp_L_ep_low_s, lfp_M_ep_low_s, lfp_H_ep_low_s, lfp_B_low_s_list, lfp_L_low_s_list, lfp_M_low_s_list, lfp_H_low_s_list)

# =============================================================================
# HIGH SPEED

# keep only good trial for high speed
lfp_B_high_s_list, lfp_L_high_s_list, lfp_M_high_s_list, lfp_H_high_s_list = keep_only_good_trials(LfpRB, LfpRL, LfpRM, LfpRH, tot_mask_B_high_s, tot_mask_L_high_s, tot_mask_M_high_s, tot_mask_H_high_s, "high speed")
# stack 1 min Lfp for high speed trials into a 4D list/array: nch, min id, id trial, length trial,
lfp_B_ep_high_s, lfp_L_ep_high_s, lfp_M_ep_high_s, lfp_H_ep_high_s = stack_lfp_1min(lfp_B_ep_high_s, lfp_L_ep_high_s, lfp_M_ep_high_s, lfp_H_ep_high_s, lfp_B_high_s_list, lfp_L_high_s_list, lfp_M_high_s_list, lfp_H_high_s_list)


print('nch ', len(lfp_B_ep_low_s), 'n min ', len(lfp_B_ep_low_s[0][0]),' size', lfp_B_ep_low_s[0][0].shape)

# =============================================================================
# Save files in matlab
# =============================================================================

save_matlab_files(rec,sess,'HPC', lfp_B_ep_low_s,lfp_L_ep_low_s,lfp_M_ep_low_s,lfp_H_ep_low_s, lfp_B_ep_high_s,lfp_L_ep_high_s,lfp_M_ep_high_s,lfp_H_ep_high_s)







