# -*- coding: utf-8 -*-
"""
Spyder Editor

@ Gino Del Ferraro, Fenton lab, Oct 2023
"""

"""
Code used to generate processed LFP for the computation of the PSD and Spectrogram
for each of the HPC subregions (e.g. Stratum Oriens, Stratum Pyramidale (Ripple band), S. Radiatum, ....).

The input signal is Neuropixel LFP (2 x N) channels in the (x, y) dimensions. 
Channels along x are averaged together after pre-processing. Channels in the y dimension
are averaged together 2-by-2. The resulting number of CSD channels is (2 x N / 4)

@ Gino Del Ferraro, Fenton lab, Oct 2023
"""
import sys
import os


base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
py_utils_path = os.path.join(base_dir, '00-functions-tools', 'python-utils-tools')

sys.path.append(py_utils_path)


from utilities_ketamine_analysis_v8 import *
from utils_signal_processing import *
from utils_plotting import *
from utils_general import *



""" Input parameters """
sess = 2 # session number 

offset = 5 # starting min for each epoch
tot_min = 30 - 2*offset # tot numb of minutes in each epoch. Each epoch starts at offset and ends at '30 min - offset'
save_var = "so_sp" # saving file name, Current Source Density

n_el_block = 1 # diving factor to get the tot numb of electrode.  n_el_block = 1 if both avg_x and avg_y are False, it's 4 if both are true, it's 2 if only one of them is True

plot_CSD = True # plot CSD and LFP for comparison
qband = 200 # Q factor in the notch filter 


""" Define Lfp channels of Subregion of HPC -- 
This is used to select only either  1. Septum Oriens; 2. Septum Oriens and Septum Pyramidale 
Either way, the selected region is called CA1 for brevity
"""
CA1_start = 0
CA1_end = 46
" ------------------- "


binFullPath = r'C:\Users\fentonlab\Desktop\Gino\LFPs'
HPC_path_file = os.path.join(r'C:\Users\fentonlab\Desktop\Gino\LFPs','HPC_lfp_paths.file')
PFC_path_file = os.path.join(r'C:\Users\fentonlab\Desktop\Gino\LFPs','PFC_lfp_paths.file')


# =============================================================================
# Load LFP and speed - upsample speed 
# =============================================================================

# ====== Load Lfp and speed data for a specific recording and brain area 
Lfp, speed, gain, rec, _, _ = load_data(binFullPath,HPC_path_file,PFC_path_file,"HPC",sess)

# ====== Detect bad (silent) Lfp channel (if it exist)
bad_flag, next_id, bad_id = detect_silent_lfp_channel(Lfp, CA1_end, 4, 4, 2500)

# ====== Upsample Speed recording
speed_up = upsample_speed(speed, Lfp, sess, LFP_rate = 2500, speed_rate = 100)

""" Get Lfp of Subregion of HPC -- This is used to select only either 
1. Septum Oriens
2. Septum Oriens and Septum Pyramidale 
Either way, the selected region is called CA1 for brevity
"""
Lfp = Lfp[:,CA1_start:CA1_end]


# =============================================================================
# PLOTTING -- sanity checks / data visualization 
# =============================================================================

# # plot speed
# plot_speed(speed_up,sess,0,100,10,fs_lfp=2500,fs_beh=100)
# # # plot speed histogram 
# plot_speed_histo_logscale(speed,sess)
# plot_speed_histo_regimes(speed_up, th_low = 30,th_mid = 100)

# # plot Lfp
# plot_lfp_two_channels(Lfp,1,2,10,100,100,N=2500)
# plot_lfp_various_channels(Lfp,0,9,0,60,3,3,10,N=2500)


# # RE-REFERENCING 
# mean_B = Lfp_B_min - np.mean(Lfp_B_min, axis=1, keepdims=True)
# plot_lfp_two_channels(Lfp_B_min,1,2,10,100,100,N=2500)
# plot_lfp_two_channels(mean_B,1,2,10,100,100,N=2500)

#%%
# ====== Split speed and Lfp into injection periods (epochs): baseline, low, mid, and high injection
Lfp_B, Lfp_L, Lfp_M, Lfp_H, speed_B, speed_L, speed_M, speed_H = split_into_epochs(Lfp,speed_up,fs=2500)
 

# nch = int(Lfp_B.shape[1]// 4) # number of channel after averaging a 2x2 block 
nch = int(Lfp_B.shape[1]// n_el_block) # number of channel after averaging a 2x2 block 

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
# all trials -- (min id, length T for 60 sec, channel id)
lfp_B_ep = []
lfp_L_ep = []
lfp_M_ep = []
lfp_H_ep = []

# mask low 
mask_B_low = []
mask_L_low = []
mask_M_low = []
mask_H_low = []

# mask high 
mask_B_high = []
mask_L_high = []
mask_M_high = []
mask_H_high = []


# =============================================================================
# SELECT ONE MINUTE DATA and iterate for 20 min
# =============================================================================


for current_min in range(0,tot_min):
    
    print('\n# ======== Current minute = {}  ----------------------- \n'.format(current_min))
    
    # =============================================================================
    # Prepare Lfp data, 1 min 
    # =============================================================================
    
    # ====== Select 1 min data 
    Lfp_B_min, Lfp_L_min, Lfp_M_min, Lfp_H_min, speed_B_min, speed_L_min, speed_M_min, speed_H_min = \
        select_1min_data(Lfp_B, Lfp_L, Lfp_M, Lfp_H, speed_B, speed_L, speed_M, speed_H, current_min, offset = 5, fs=2500)
    
    # ====== Replace Lfp bad channel with nearest neighbor (if bad channel exists)
    if bad_flag:
        Lfp_B_min, Lfp_L_min, Lfp_M_min, Lfp_H_min = replace_bad_lfp_channel(Lfp_B_min, Lfp_L_min, Lfp_M_min, Lfp_H_min, bad_id, next_id)
   
    # === plot bad channel replaced and nearest neighbor
    # plot_lfp_two_channels(Lfp_L_min,bad_id,next_id,0,60,10,N=2500)

    # ====== Average Lfp in Neuropixel at the same depth (avg 2 electrodes together)
       
    
    # ====== Average Lfp in Neuropixelin a 2x2 channel block (avg 4 electrodes together)
    # Lfp_B_avg, Lfp_L_avg, Lfp_M_avg, Lfp_H_avg = average_lfp_4_channels(Lfp_B_min,Lfp_L_min,Lfp_M_min,Lfp_H_min)
    # ====== Do not perform any average, this is to be used in single cell analysis 
    Lfp_B_avg, Lfp_L_avg, Lfp_M_avg, Lfp_H_avg = Lfp_B_min, Lfp_L_min, Lfp_M_min, Lfp_H_min 


    # =============================================================================
    # Filter 1 min LFP (band pass)
    # =============================================================================
    
    # ====== Filter Lfp in each epoch
    lfp_filt_B, lfp_filt_L, lfp_filt_M, lfp_filt_H = filter_lfp_in_each_epoch(Lfp_B_avg, Lfp_L_avg, Lfp_M_avg, Lfp_H_avg, gain, qband)
          
    # ====== Decimate Lfp and speed (subsample): input at 2500 Hz, output at 1250 Hz
    
    # ====== Decimate Lfp and speed (subsample)
    lfp_dec_B, lfp_dec_L, lfp_dec_M, lfp_dec_H, speed_dec_B, speed_dec_L, speed_dec_M, speed_dec_H = \
        decimate_lfp_and_speed(lfp_filt_B, lfp_filt_L, lfp_filt_M, lfp_filt_H, speed_B_min,speed_L_min,speed_M_min,speed_H_min)
    
    # ====== Stack lfp all trials for each minute together, (min id, length T for 60 sec, channel id) 
    lfp_B_ep,lfp_L_ep,lfp_M_ep,lfp_H_ep = \
        stack_lfp_1min_all_trials(lfp_B_ep, lfp_L_ep, lfp_M_ep, lfp_H_ep, lfp_dec_B, lfp_dec_L, lfp_dec_M, lfp_dec_H)



    # =============================================================================
    # MASKING SPEED AND LFP ARTIFACTS 
    # =============================================================================
    
    # ====== Combine speed mask with LFP artifacts mask 
    tot_mask_B_low_s, tot_mask_L_low_s, tot_mask_M_low_s, tot_mask_H_low_s, tot_mask_B_high_s, tot_mask_L_high_s, tot_mask_M_high_s,tot_mask_H_high_s = \
        make_speed_and_lfp_maks(lfp_dec_B,lfp_dec_L, lfp_dec_M, lfp_dec_H, speed_dec_B, speed_dec_L, speed_dec_M, speed_dec_H, win = 1250, th = 30)
  
    # ====== Stack masks all trials for each minute together
    mask_B_low, mask_L_low, mask_M_low, mask_H_low, mask_B_high, mask_L_high, mask_M_high, mask_H_high = \
        stack_mask_1min(mask_B_low, mask_L_low, mask_M_low, mask_H_low, 
                         mask_B_high, mask_L_high, mask_M_high, mask_H_high,
                         tot_mask_B_low_s, tot_mask_L_low_s, tot_mask_M_low_s, tot_mask_H_low_s,
                         tot_mask_B_high_s, tot_mask_L_high_s, tot_mask_M_high_s, tot_mask_H_high_s)
    
    
    # =============================================================================
    #  Reshape Lfp into: trial number, trial length, channels
    # =============================================================================
    
    # ====== reshape Lfp in trial x time window x channels
    win = 1250
    LfpRB = lfp_dec_B.reshape(int(lfp_dec_B.shape[0]/win),-1,lfp_dec_B.shape[1]) # reshape Lfp: trial x win length x channel
    LfpRL = lfp_dec_L.reshape(int(lfp_dec_L.shape[0]/win),-1,lfp_dec_L.shape[1]) # reshape Lfp: trial x win length x channel
    LfpRM = lfp_dec_M.reshape(int(lfp_dec_M.shape[0]/win),-1,lfp_dec_M.shape[1]) # reshape Lfp: trial x win length x channel
    LfpRH = lfp_dec_H.reshape(int(lfp_dec_H.shape[0]/win),-1,lfp_dec_H.shape[1]) # reshape Lfp: trial x win length x 
    
    
    
    print('reshaped Lfp, ', LfpRB.shape, LfpRL.shape, LfpRM.shape, LfpRH.shape)
    
    
    # =============================================================================
    # Keep good trials only and stack them into a 4D array
    # =============================================================================
    
    # =============================================================================
    # LOW SPEED trials 

    # ====== keep only good trial for low speed
    lfp_B_low_s_list, lfp_L_low_s_list, lfp_M_low_s_list, lfp_H_low_s_list = \
        keep_only_good_trials(LfpRB, LfpRL, LfpRM, LfpRH, tot_mask_B_low_s, tot_mask_L_low_s, tot_mask_M_low_s, tot_mask_H_low_s, "low speed")
    # ====== stack 1 min Lfp for low speed trials into a 4D list/array: nch, min id, id trial, length trial,
    lfp_B_ep_low_s, lfp_L_ep_low_s, lfp_M_ep_low_s, lfp_H_ep_low_s = \
        stack_lfp_1min(lfp_B_ep_low_s, lfp_L_ep_low_s, lfp_M_ep_low_s, lfp_H_ep_low_s, lfp_B_low_s_list, lfp_L_low_s_list, lfp_M_low_s_list, lfp_H_low_s_list)
    
    # =============================================================================
    # HIGH SPEED trials 
    
    # ====== keep only good trial for high speed
    lfp_B_high_s_list, lfp_L_high_s_list, lfp_M_high_s_list, lfp_H_high_s_list = \
        keep_only_good_trials(LfpRB, LfpRL, LfpRM, LfpRH, tot_mask_B_high_s, tot_mask_L_high_s, tot_mask_M_high_s, tot_mask_H_high_s, "high speed")
    # ====== stack 1 min Lfp for high speed trials into a 4D list/array: nch, min id, id trial, length trial,
    lfp_B_ep_high_s, lfp_L_ep_high_s, lfp_M_ep_high_s, lfp_H_ep_high_s = \
        stack_lfp_1min(lfp_B_ep_high_s, lfp_L_ep_high_s, lfp_M_ep_high_s, lfp_H_ep_high_s, lfp_B_high_s_list, lfp_L_high_s_list, lfp_M_high_s_list, lfp_H_high_s_list)
    
    


# =============================================================================
# Save files in matlab
# # low speed and high speed trials 
# print('Saving Lfp split into trials ...')
# save_matlab_files(rec, sess, 'HPC', 
#                   lfp_B_ep_low_s, lfp_L_ep_low_s, lfp_M_ep_low_s, lfp_H_ep_low_s, 
#                   lfp_B_ep_high_s, lfp_L_ep_high_s, lfp_M_ep_high_s, lfp_H_ep_high_s)


# All trials with mask for low/high speed --- ONLY FOR CA1 -- For single LFP-spike cell analysis 
# (min id, length T for 60 sec, channel id)



"""  Save LFP trials without artifacts for low and high speed - usage: PSD calculation """

print('Saving Lfp split into trials ...')
save_matlab_files(rec, sess, 'HPC', 
                  lfp_B_ep_low_s, lfp_L_ep_low_s, lfp_M_ep_low_s, lfp_H_ep_low_s, 
                  lfp_B_ep_high_s, lfp_L_ep_high_s, lfp_M_ep_high_s, lfp_H_ep_high_s, save_var)



# data_B_low = load_lfp_data(r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC\2022-08-01_04-30-00_M015_RSK_mPFC_HPC_3_10_30mpk\lfp_B_epoch_low_speed.mat')
# data_L_low = load_lfp_data(r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC\2022-08-01_04-30-00_M015_RSK_mPFC_HPC_3_10_30mpk\lfp_L_epoch_low_speed.mat')

""" Save all LFP trials and total mask - usage: Spectrograms, and other signal processing analysis """

print('Saving lfp whole min recording + masks ...')
save_matlab_files_all_lfps(rec,sess,'HPC', lfp_B_ep, lfp_L_ep, lfp_M_ep, lfp_H_ep, 
                               mask_B_low, mask_L_low, mask_M_low, mask_H_low, 
                               mask_B_high, mask_L_high, mask_M_high, mask_H_high, save_var)



