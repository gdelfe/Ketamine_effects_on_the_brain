# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 16:36:43 2023

@author: Gino Del Ferraro, Fenton Lab, Oct 2023
"""

import numpy as np
import pickle
from DemoReadSGLXData.readSGLX import readMeta, SampRate, makeMemMapRaw, ExtractDigital, Int2Volts
from pathlib import Path

# import ghostipy as gsp
import numpy as np
import scipy.signal as signal
from scipy.io import savemat
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import pandas as pd
from IPython.display import display
import os
import h5py
import copy
from scipy.stats import zscore
from urllib.request import urlopen
plt.rcParams['pdf.fonttype'] = 42

#import mne
import sys
import pickle
from utils_plotting import *


def load_data(binFullPath,HPC_path_file,PFC_path_file,brain_reg,sess):
    
    
    rec = [{"PFC":"path","HPC":"path"} for _ in range(36)] # initialize dictonary for all the paths
    
    # =============================================================================
    #     LOAD file paths 
    # =============================================================================
    # if brain region is HPC load HPC paths 
    if brain_reg == 'HPC':
    # load HPC file names and store them in rec
        with open(HPC_path_file,'rb') as f:
            HPC_file_list = pickle.load(f)
            
        for isess, path in enumerate(HPC_file_list):
            rec[isess]['HPC'] =  str(path)
    
    elif brain_reg == 'PFC':
        # load PFC file names and store them in rec 
        with open(PFC_path_file,'rb') as f:
            PFC_file_list = pickle.load(f)
            
        for isess, path in enumerate(PFC_file_list):
            rec[isess]['PFC'] =  str(path)
        
    else:
        sys.exit('Brain region is not HPC or PFC -- exit')
        
        
    # select path for specific session recording 
    binFullPath = rec[sess][brain_reg] 
    print('Loading file in: ',binFullPath)
    # print(binFullPath)
    numChannels = 385
    
    # =============================================================================
    # # load  LFP data 
    # =============================================================================
    
    rawDataLfp = np.memmap(binFullPath, dtype='int16', mode='r')
    dataLfp = np.reshape(rawDataLfp, (int(rawDataLfp.size/numChannels), numChannels))
    
    ### Load scaling factor for Lfp
    meta = readMeta(Path(binFullPath))
    
    # 3A, 3B1, 3B2 (NP 1.0)
    imroList = meta['imroTbl'].split(sep=')')[1:-1]
    LFgain = np.array([int(i.split(sep=' ')[4]) for i in imroList])
    fI2V = Int2Volts(meta)
    conv = fI2V/LFgain
    gain = np.append(conv, 1) # scaling factor for Lfp
    
    ### Check if gain_correct is always the same 
    if np.unique(gain[0:-1]):
        print(f'Gain factors all the same for each Lfp channel: {gain[0]} in Volts')
        gain = gain[0] 
    else:
        sys.exit('Gain factor multiplication should be done channel by channel, this is just a bit slower')
        
        
    #### Read LFP start/end time, HPC channels - Choose whether loading HPC or PFC here
        
    pd.set_option('display.max_colwidth',100)
    in_file = os.path.join(r'C:\Users\fentonlab\Desktop\Gino\LFPs\CA1_DG_id', "recid_CA1_DG_id.csv")
    Lfp_aln = pd.read_csv(in_file)
    

    print(Lfp_aln.head())
    
    
    # =============================================================================
    #  Trim LFP 
    # =============================================================================
    #### Pull Lfp 'start' and 'end' for the specific recording 
    
    start = Lfp_aln["Lfp start"][sess] # starting of the trial in Lfp data points
    end = Lfp_aln["Lfp end"][sess] # end of trial
    
    # multiply by two to account for the full channel list 
    ch_start = Lfp_aln["theta_start"][sess]*2  # Hpc first channel
    ch_end = Lfp_aln["theta_end"][sess]*2  # Hpc last channel
    

    # Trim Lfp to specific channels in the HPC and to specific start/end times 
    Lfp= dataLfp[start:end,ch_start:ch_end]
    # Lfp = dataLfp[start:end,:]
    print('Lfp shape after trimming ',Lfp.shape)
        
    # =============================================================================
    # LOAD SPEED 
    # =============================================================================
    
    ### Load Speed, x, and y
    speed_path = r'C:\Users\fentonlab\Desktop\Gino\speed\\'
    # x = np.load(os.path.join(speed_path, "x_aln.npy"), allow_pickle=True)
    # y = np.load(os.path.join(speed_path, "y_aln.npy"), allow_pickle=True)
    speed = np.load(os.path.join(speed_path, "speed_aln.npy"), allow_pickle=True)
    
    #### difference in time alignment between behavior and Lfp
    print('Length in minutes of speed recordings: ',(speed[sess].size)/100/60)
    print('Length in minutes of Lfp recordings: ',(Lfp.shape[0])/2500/60)
    print('Difference in seconds: ',((Lfp.shape[0])/2500/60 - (speed[sess].size)/100/60)*60)

    return Lfp, speed, gain, rec



# =============================================================================


def detect_silent_lfp_channel(Lfp, N = 2500):
    
    current_min = 0 # current min to look at 
    length = 3 # length of period to look at, i.e. 1 = 1 min 
    offset = 5 # starting of epoch in min

    start = (offset + current_min)*60*N - 1    # start point in time points
    end = (offset + current_min + length)*60*N - 1   # stop point in time points

    # Lfp trim: 3 minute, one every M channel
    Lfp_B_3min = Lfp[start:end,:]  # base line period
    

    # remove mean from Lfp for each channel
    lfp_ms = Lfp[start:end,:] - np.mean(Lfp[start:end,:],axis=0)
    # take abs value for each lfp channel
    lfp_abs_all = np.abs(lfp_ms)
    # compute mean across time for abs lfp
    mean_abs = np.mean(lfp_abs_all,axis=0)
    
    min_val = np.min(mean_abs)
    min_id = np.argmin(mean_abs)
    
    print('Bad channel Lfp abs average over 3 min time: ',min_val, 'Bad channel ID: ', min_id)
    
    
    if min_id % 2: # if odd channel 
        next_id = min_id - 1
    else:           # if even channel
        next_id = min_id + 1 
    
    # plot bad channel and channel next to it
    plot_lfp_two_channels_together(Lfp,next_id,min_id,10,200,20,N=2500)
    
    return 
    

# =============================================================================

# Average 2 channels in the Neuropixel array which are on the same depth, i.e.
# average consecutive channels in the Lfp map, for each epoch separately, for 1 min of data at the time 

def average_lfp_same_depth(Lfp_B_min,Lfp_L_min,Lfp_M_min,Lfp_H_min):
    
    print("Averaging lfp ...")
    
    if int(Lfp_B_min.shape[1]/2) % 2:
        sys.exit("Number of channels in the brain region considered is not EVEN! Check channel list!")
        
    # baseline
    Lfp_RS_B = Lfp_B_min.reshape(Lfp_B_min.shape[0],int(Lfp_B_min.shape[1]/2),2) # reshape Lfp, such that channels on the same depth are into adjacent columns
    Lfp_avg_B = Lfp_RS_B.mean(axis=2)
    # low injection 
    Lfp_RS_L = Lfp_L_min.reshape(Lfp_L_min.shape[0],int(Lfp_L_min.shape[1]/2),2) # reshape Lfp, such that channels on the same depth are into adjacent columns
    Lfp_avg_L = Lfp_RS_L.mean(axis=2)
    # mid injection 
    Lfp_RS_M = Lfp_M_min.reshape(Lfp_M_min.shape[0],int(Lfp_M_min.shape[1]/2),2) # reshape Lfp, such that channels on the same depth are into adjacent columns
    Lfp_avg_M = Lfp_RS_M.mean(axis=2)
    # high injection 
    Lfp_RS_H = Lfp_B_min.reshape(Lfp_H_min.shape[0],int(Lfp_H_min.shape[1]/2),2) # reshape Lfp, such that channels on the same depth are into adjacent columns
    Lfp_avg_H = Lfp_RS_H.mean(axis=2)
    
    return Lfp_avg_B, Lfp_avg_L, Lfp_avg_M, Lfp_avg_H

# =============================================================================

def split_into_epochs(Lfp,speed_up,N=2500):
    
    # cut Lfp and speed up to 2 h time window (disregard data above 2 h)
    L = N*60*30*4 # time points in a 2 h time window
    # trim Lfp and speed 
    speed = speed_up[0:L]
    lfp = Lfp[0:L,:]
    
    speed_periods = speed.reshape(-1,int(speed.size/4))  # reshape speed in baseline, low, mid, high injection time periods
    
    win_30 = N*60*30 # 30 min window 
    Lfp_B = lfp[0:win_30,:]  # base line period
    Lfp_L = lfp[win_30:2*win_30,:]  # low injection 
    Lfp_M = lfp[2*win_30:3*win_30,:]  # mid injection 
    Lfp_H = lfp[3*win_30:4*win_30,:]  # high injection 
    
    speed_B = speed_periods[0,:]
    speed_L = speed_periods[1,:]
    speed_M = speed_periods[2,:]
    speed_H = speed_periods[3,:]
    
    print('min in each epoch: ',speed_B.size/N/60)
    
    return Lfp_B, Lfp_L, Lfp_M, Lfp_H, speed_B, speed_L, speed_M, speed_H

# =============================================================================


def select_1min_data(Lfp_B, Lfp_L, Lfp_M, Lfp_H, speed_B, speed_L, speed_M, speed_H, current_min, N=2500):
    
    #### 3. Select 1 min of data for each Epoch, for both Lfp and speed
    ##### Select 1 every M channels for the Lfp


    M = 1 # keep every M channel
    length = 1 # length of period to look at, i.e. 1 = 1 min 
    offset = 5 # starting of epoch in min

    start = (offset + current_min)*60*N - 1    # start point in time points
    end = (offset + current_min + length)*60*N - 1   # stop point in time points

    # Lfp trim: 1 minute, one every M channel
    Lfp_B_min = Lfp_B[start:end,::M]  # base line period
    Lfp_L_min = Lfp_L[start:end,::M]  # low injection 
    Lfp_M_min = Lfp_M[start:end,::M]  # mid injection 
    Lfp_H_min = Lfp_H[start:end,::M]  # high injection 

    # speed
    speed_B_min = speed_B[start:end] 
    speed_L_min = speed_L[start:end]
    speed_M_min = speed_M[start:end]
    speed_H_min = speed_H[start:end]

    print('1 min data: Lfp shape {}, speed shape {}, length in sec: {}\n'.format(Lfp_B_min.shape, speed_B_min.shape, speed_B_min.size/N))
    
    return Lfp_B_min, Lfp_L_min, Lfp_M_min, Lfp_H_min, speed_B_min, speed_L_min, speed_M_min, speed_H_min
         

# =============================================================================

    
def decimate_lfp_and_speed(lfp_filt_B,lfp_filt_L,lfp_filt_M,lfp_filt_H,speed_B_min,speed_L_min,speed_M_min,speed_H_min):
    
    
    # decimated Lfp to 1250 Hz
    lfp_dec_B = lfp_filt_B[::2,:]
    lfp_dec_L = lfp_filt_L[::2,:]
    lfp_dec_M = lfp_filt_M[::2,:]
    lfp_dec_H = lfp_filt_H[::2,:]
    
    # decimated speed to 1250 Hz
    speed_dec_B = speed_B_min[::2]
    speed_dec_L = speed_L_min[::2]
    speed_dec_M = speed_M_min[::2]
    speed_dec_H = speed_H_min[::2]
    
    return lfp_dec_B, lfp_dec_L, lfp_dec_M, lfp_dec_H, speed_dec_B, speed_dec_L, speed_dec_M, speed_dec_H

# =============================================================================


# input: 
# lfp_dec: lfp (subsampled)
# win: time length of the window we are looking at in time points, i.e. 1 sec = 1250 points
# std_th: std threshold
# period: string for baseline, low, etc..

def lfp_artifacts_mask(lfp_dec,win,std_th):

    zlfp = zscore(lfp_dec,axis=0)
    mask_lfp = np.abs(zlfp) < std_th

    nch = lfp_dec.shape[1] # num of channels
    T_length = lfp_dec.shape[0] # time series length 

    mask_trial = [[] for ch in range(nch)] # create mask for each channel, per trial
    mask = [[] for ch in range(nch)] # create mask for each channel, per time point  
    good_trial_cnt = []
    
    for ch in range(nch): # for each channel 
        cnt = 0 # count 
        for i in range(0,T_length,win): # for time length
            data_win = mask_lfp[i:i+win,ch]
            if np.all(data_win): # if there are no artifact within the win considered
                mask_trial[ch].extend([True]) # create True mask for that window, per trial
                mask[ch].extend([1]*win) # create True mask for that window, per time point 
                cnt +=1
            else:
                mask_trial[ch].extend([False]) # create False mask for that window, per trial
                mask[ch].extend([-1]*win) # create True mask for that window, per time point 
        
        good_trial_cnt.append(cnt/(T_length/win)) # cnt n. of good trials in each channel 

    mask_trial_arr = np.array(mask_trial).T
    mask_arr = np.array(mask).T
    good_trial_arr = np.array(good_trial_cnt)
    
    return mask_trial_arr, mask_arr, good_trial_arr


# =============================================================================

# Create a speed mask. 

# True if speed < threshold at any point in a given time window of length win
# False if there is at least one value above threshold in the given time window
# Input: 
# speed_up: upsampled speed
# win: amplitude in data points of the win we are looking at for thresholding
# thresh: speed threshold for speed
# period: string for baseline, low, mid, high injection 

def create_speed_mask(speed_up,win,thresh,level,period):
    mask = []
    if level == 'low':
        speed_th = speed_up < thresh # speed mask below threshold
    elif level == 'high':
        speed_th = speed_up > thresh # speed mask below threshold
    else:
        sys.exit('Wrong speed scenario, choose either low or high for level')
            
    cnt = 0 # count for low speed trials 
    for i in range(0,len(speed_th),win):
        data_win = speed_th[i:i+win]
        if np.all(data_win): # if all the speed value are below threshold
            mask.extend([True]) # create True mask for that window
            cnt +=1
        else:
            mask.extend([False]) # create False mask for that window
    print("{} speed trials % in {} : {:.2f}".format(level,period,cnt/((len(speed_th)/win))))
    return np.array(mask)


# =============================================================================


# Create mask for:
#     1. Lfp artifacts
#     2. speed (low and high speed)
#     3. Combine them together

# win = window length, e.g. 1250 = 1 sec
# th = speed threshold 
    
def make_speed_and_lfp_maks(lfp_dec_B,lfp_dec_L, lfp_dec_M, lfp_dec_H, speed_dec_B, speed_dec_L, speed_dec_M, speed_dec_H, win = 1250, th = 30):
    
    nch = lfp_dec_B.shape[1] # numb of channels 
    
    # =============================================================================
    #     LFP mask for artifacts 
    # =============================================================================
        
    # Create mask for artifact trials in Lfp

    mask_trial_B, lfp_mask_B, good_trial_rate_B = lfp_artifacts_mask(lfp_dec_B,win,4)
    mask_trial_L, lfp_mask_L, good_trial_rate_L = lfp_artifacts_mask(lfp_dec_L,win,4)
    mask_trial_M, lfp_mask_M, good_trial_rate_M = lfp_artifacts_mask(lfp_dec_M,win,4)
    mask_trial_H, lfp_mask_H, good_trial_rate_H = lfp_artifacts_mask(lfp_dec_H,win,4)

    # =============================================================================
    # Speed masks for Low and High speed
    # =============================================================================

    # make a speed map for low speed
    speed_mask_B_low = create_speed_mask(speed_dec_B,win,th,'low','base')
    speed_mask_L_low = create_speed_mask(speed_dec_L,win,th,'low','low')
    speed_mask_M_low = create_speed_mask(speed_dec_M,win,th,'low','mid')
    speed_mask_H_low = create_speed_mask(speed_dec_H,win,th,'low','high')
    print()
    # make a speed map for low speed
    speed_mask_B_high = create_speed_mask(speed_dec_B,win,th,'high','base')
    speed_mask_L_high = create_speed_mask(speed_dec_L,win,th,'high','low')
    speed_mask_M_high = create_speed_mask(speed_dec_M,win,th,'high','mid')
    speed_mask_H_high = create_speed_mask(speed_dec_H,win,th,'high','high')
    print()
    
    # =============================================================================
    # Combine speed mask and Lfp artifact mask
    # =============================================================================
    
    # Combine the lfp good trial mask with the low speed mask
    tot_mask_B_low_s = (mask_trial_B.T & speed_mask_B_low).T
    tot_mask_L_low_s = (mask_trial_L.T & speed_mask_L_low).T
    tot_mask_M_low_s = (mask_trial_M.T & speed_mask_M_low).T
    tot_mask_H_low_s = (mask_trial_H.T & speed_mask_H_low).T
    
    # cnt = np.sum(mask_trial_L,axis=0)
    # cnt_tot = np.sum(tot_mask_L_low_s,axis=0)
    # cnt - cnt_tot # mask difference, must be >= 0
    
    # print('\nNo artifact trials\n ',cnt)
    # print('\nNo artifact, low speed trials \n ',cnt_tot)
    
    # Combine the lfp good trial mask with the high speed mask
    tot_mask_B_high_s = (mask_trial_B.T & speed_mask_B_high).T
    tot_mask_L_high_s = (mask_trial_L.T & speed_mask_L_high).T
    tot_mask_M_high_s = (mask_trial_M.T & speed_mask_M_high).T
    tot_mask_H_high_s = (mask_trial_H.T & speed_mask_H_high).T
    
    # cnt = np.sum(mask_trial_M,axis=0)
    # cnt_tot = np.sum(tot_mask_M_high,axis=0)
    # cnt - cnt_tot # mask difference, must be >= 0
    
    # print('\nNo artifact trials\n ',cnt)
    # print('\nNo artifact, high speed trials \n ',cnt_tot)
    
    # low speed trial rate
    tot_good_trial_rate_B = np.sum(np.sum(tot_mask_B_low_s,axis=0))/60/nch
    tot_good_trial_rate_L = np.sum(np.sum(tot_mask_L_low_s,axis=0))/60/nch
    tot_good_trial_rate_M = np.sum(np.sum(tot_mask_M_low_s,axis=0))/60/nch
    tot_good_trial_rate_H = np.sum(np.sum(tot_mask_H_low_s,axis=0))/60/nch
    print('good trial-low speed rate, base = {:.2f}, low = {:.2f}, mid = {:.2f}., high = {:.2f}'.format(tot_good_trial_rate_B,tot_good_trial_rate_L,tot_good_trial_rate_M,tot_good_trial_rate_H))
    
    # high speed trial rate
    tot_good_trial_rate_B = np.sum(np.sum(tot_mask_B_high_s,axis=0))/60/nch
    tot_good_trial_rate_L = np.sum(np.sum(tot_mask_L_high_s,axis=0))/60/nch
    tot_good_trial_rate_M = np.sum(np.sum(tot_mask_M_high_s,axis=0))/60/nch
    tot_good_trial_rate_H = np.sum(np.sum(tot_mask_H_high_s,axis=0))/60/nch
    print('good trial-high speed rate, base = {:.2f}, low = {:.2f}, mid = {:.2f}., high = {:.2f}'.format(tot_good_trial_rate_B,tot_good_trial_rate_L,tot_good_trial_rate_M,tot_good_trial_rate_H))
    print()
    
    return tot_mask_B_low_s, tot_mask_L_low_s, tot_mask_M_low_s, tot_mask_H_low_s, tot_mask_B_high_s, tot_mask_L_high_s, tot_mask_M_high_s,tot_mask_H_high_s 


# =============================================================================

def keep_only_good_trials(LfpRB,LfpRL,LfpRM,LfpRH, tot_mask_B,tot_mask_L,tot_mask_M,tot_mask_H, speed_string):
    
    
    # Keep only good trial in Lfp for 1 min recording
    # ouput: list with nch, each element contains good-trial-idx x time values, good-trial-idx may differ channel by channel
    
    lfp_B_list = [] 
    lfp_L_list = []
    lfp_M_list = []
    lfp_H_list = []  
    
    nch = LfpRB.shape[2] # numb of channels 
    for ch in range(nch):
        good_trials = LfpRB[tot_mask_B[:,ch],:,ch] # good trials per each channel
        lfp_B_list.append(good_trials)
        good_trials = LfpRL[tot_mask_L[:,ch],:,ch] # good trials per each channel
        lfp_L_list.append(good_trials)
        good_trials = LfpRM[tot_mask_M[:,ch],:,ch] # good trials per each channel
        lfp_M_list.append(good_trials)
        good_trials = LfpRH[tot_mask_H[:,ch],:,ch] # good trials per each channel
        lfp_H_list.append(good_trials)
    
    print('\n',speed_string,', nch:', len(lfp_B_list), lfp_B_list[0].shape, lfp_L_list[0].shape, lfp_L_list[0].shape, lfp_L_list[0].shape)
        
        
    return lfp_B_list, lfp_L_list, lfp_M_list, lfp_H_list

# =============================================================================

def stack_lfp_1min(lfp_B_epoch,lfp_L_epoch,lfp_M_epoch,lfp_H_epoch,lfp_B_list,lfp_L_list,lfp_M_list,lfp_H_list):
    
    nch = len(lfp_B_list)
    for ch in range(nch):
        lfp_B_epoch[ch].append(lfp_B_list[ch]) # add one minute lfp: trial num x trial length, for each channel
        lfp_L_epoch[ch].append(lfp_L_list[ch])
        lfp_M_epoch[ch].append(lfp_M_list[ch])
        lfp_H_epoch[ch].append(lfp_H_list[ch])
    

    return lfp_B_epoch, lfp_L_epoch, lfp_M_epoch, lfp_H_epoch

# =============================================================================


def save_matlab_files(rec,sess,brain_reg, lfp_B_ep_low_s,lfp_L_ep_low_s,lfp_M_ep_low_s,lfp_H_ep_low_s, lfp_B_ep_high_s,lfp_L_ep_high_s,lfp_M_ep_high_s,lfp_H_ep_high_s):
    
    main_dir = r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC'
    path = rec[sess][brain_reg]
    
    dir_sess = path.split('\\')[-3] # path for session directory
    full_dir_path = os.path.join(main_dir,dir_sess)
    
    if not os.path.exists(full_dir_path):
        os.makedirs(full_dir_path)
    
    # =============================================================================
    #     Low Speed Lfp 
    # =============================================================================
    
    # file path to save in matlab 
    out_file_B = os.path.join(full_dir_path, "lfp_B_epoch_low_speed.mat")
    out_file_L = os.path.join(full_dir_path, "lfp_L_epoch_low_speed.mat")
    out_file_M = os.path.join(full_dir_path, "lfp_M_epoch_low_speed.mat")
    out_file_H = os.path.join(full_dir_path, "lfp_H_epoch_low_speed.mat")
    
    # create dictionaries to save in matlab
    mat_lfp_B = {'lfp_B': lfp_B_ep_low_s}
    mat_lfp_L = {'lfp_L': lfp_B_ep_low_s}
    mat_lfp_M = {'lfp_M': lfp_B_ep_low_s}
    mat_lfp_H = {'lfp_H': lfp_B_ep_low_s}
    
    # save lfp for each epoch in matlab format 
    savemat(out_file_B,mat_lfp_B)
    savemat(out_file_L,mat_lfp_L)
    savemat(out_file_M,mat_lfp_M)
    savemat(out_file_H,mat_lfp_H)
    
    # =============================================================================
    #     High Speed Lfp 
    # =========================================================================
    
    # file path to save in matlab 
    out_file_B = os.path.join(full_dir_path, "lfp_B_epoch_high_speed.mat")
    out_file_L = os.path.join(full_dir_path, "lfp_L_epoch_high_speed.mat")
    out_file_M = os.path.join(full_dir_path, "lfp_M_epoch_high_speed.mat")
    out_file_H = os.path.join(full_dir_path, "lfp_H_epoch_high_speed.mat")
    
    # create dictionaries to save in matlab
    mat_lfp_B = {'lfp_B': lfp_B_ep_high_s}
    mat_lfp_L = {'lfp_L': lfp_B_ep_high_s}
    mat_lfp_M = {'lfp_M': lfp_B_ep_high_s}
    mat_lfp_H = {'lfp_H': lfp_B_ep_high_s}
    
    # save lfp for each epoch in matlab format 
    savemat(out_file_B,mat_lfp_B)
    savemat(out_file_L,mat_lfp_L)
    savemat(out_file_M,mat_lfp_M)
    savemat(out_file_H,mat_lfp_H)
    
    
    return 





