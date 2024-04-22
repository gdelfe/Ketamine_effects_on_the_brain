# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 16:36:43 2023

@author: Gino Del Ferraro, Fenton Lab, Oct 2023
"""
from utilities_ketamine_analysis_v8 import *

from scipy.io import savemat
import matplotlib.pyplot as plt
import os
from scipy.stats import zscore
plt.rcParams['pdf.fonttype'] = 42

import sys
from utils_plotting import *
import h5py


def read_meta(binpath):
    metapath = binpath[:-len('.bin')] + ".meta"
    md = {}
    assert os.path.isfile(metapath)
    with open(metapath, 'r') as f:
        lines = f.read().splitlines()
        # convert the list entries into key value pairs
        for m in lines:
            vals = m.split(sep='=')
            if vals[0][0] == '~':
                k = vals[0][1:len(vals[0])]
            else:
                k = vals[0]
            md.update({k: vals[1]})
    return(md)


# Return a multiplicative factor for converting 16-bit file data
# to voltage. This does not take gain into account. The full
# conversion with gain is:
#         dataVolts = dataInt * fI2V / gain
# Note that each channel may have its own gain.
def get_i2v(meta):
    if meta['typeThis'] == 'imec':
        if 'imMaxInt' in meta:
            maxint = int(meta['imMaxInt'])
        else:
            maxint = 512
        fI2V = float(meta['imAiRangeMax']) / maxint
    else:
        fI2V = float(meta['niAiRangeMax']) / 32768
    return(fI2V)


"""
Load LFP data, trim it in order to align it with speed data
Load speed data 
Output: 
    Lfp, speed, gain (Lfp scaling factor), rec (list of session names), 
    ch_start (start channel for HPC), ch_end (end channel for HPC)
"""
def load_session():
    # select path for specific session recording
    root = '/Users/lukearend/phd/kfx/data'
    file = '2022-08-11-01-55-00_M018_SAL_mPFC_HPC_0_0_0mpk_g0_t0.imec1.lf.bin'
    speedpath = '/Users/lukearend/phd/kfx/data/behaviour/speed_aln.npy'
    hpc_start = 60 * 2
    hpc_end = 123 * 2
    lfp_start = 38439
    lfp_end = 18060209

    nchan = 385

    path = os.path.join(root, file)
    meta = read_meta(path)
    speed = np.load(speedpath, allow_pickle=True)[7]

    print('Loading ', path)
    lfp = np.memmap(path, dtype='int16', mode='r')
    lfp = np.reshape(lfp, (int(lfp.size / nchan), nchan))

    # 3A, 3B1, 3B2 (NP 1.0)
    gaintable = meta['imroTbl'].split(sep=')')[1:-1]
    gain = np.array([int(i.split(sep=' ')[4]) for i in gaintable])
    assert np.all(gain[0:-1] == gain[0])
    i2v = float(meta['imAiRangeMax']) / int(meta['imMaxInt'])
    i2v = i2v / gain[0] # int16 to Volts conversion factor


    # Trim LFP to specific channels in the HPC and to specific start/end times
    lfp = lfp[lfp_start:lfp_end,hpc_start:hpc_end]
    print('LFP shape after trimming: ', lfp.shape)

    # difference in time alignment between behavior and LFP
    print('Length in minutes of speed recordings: ', (speed.shape[0]) / 100 / 60)
    print('Length in minutes of LFP recordings: ', (lfp.shape[0]) / 2500 / 60)
    print('Difference in seconds: ', speed.shape[0] / 100 - lfp.shape[0] / 2500)
    return speed, lfp, i2v


"""
Detect bad (silent) Lfp channel. 

The bad Lfp channel has much less activity than 
any other channel (about 1 order of magnitude less)
Procedure:
Scale each lfp channel by its mean, take the abs for each time point, take the min
of the abs(lfp), find the channel with the minimum abs(lfp) and flag it as bad channel
"""
def detect_bad_channel(lfp, threshold=4):
    start = 5 * 60 * 2500  # start point in time points
    end = 8 * 60 * 2500 # stop point in time points

    # remove mean from lfp for each channel
    lfp_ms = lfp[start:end,:] - lfp[start:end,:].mean(axis=0)
    # take abs value for each lfp channel
    lfp_abs_all = np.abs(lfp_ms)
    # compute mean across time for abs lfp
    mean_abs = lfp_abs_all.mean(axis=0)
    mean_tot = mean_abs.mean()
    print(f'Total mean abs LFP across channels over time 05:00-08:00: {mean_tot:.4f}')

    min_val = np.min(mean_abs)
    bad_id = np.argmin(mean_abs)

    # if bad channel detected
    if min_val < mean_tot / threshold:
        next_id = bad_id - 1 if bad_id % 2 else bad_id + 1
        print(f'Bad channel LFP abs average over time 05:00-08:00: {min_val:.2f}')
        print('Bad channel ID:', bad_id)
        print('Nearest neighbor channel to bad channel:', next_id)
        return bad_id, next_id

    # there is no bad channel
    print('No bad channel detected\n')
    return None, None


""" Upsample speed data to 2500 Hz, i.e. same resolution as LFP data """
def upsample_speed(speed, lfp):
    # upsample speed to match LFP, take first two hours
    t1 = np.arange(0, len(speed)) / 100
    f1 = interp1d(t1, speed, kind='linear', fill_value="extrapolate")
    t2 = np.arange(0, len(lfp)) / 2500
    speed = f1(t2)
    print(f'Upsampled speed shape {speed.shape}, LFP shape {lfp.shape}')
    return speed


""" Split Lfp and speed into 30 min Epochs: baseline, low dose injection, mid dose, high dose """
def split_into_epochs(lfp, speed, fs=2500):
    # disregard data above 2 h
    n = 4 * 30 * 60 * fs
    speed = speed[0:n]
    lfp = lfp[0:n, :]

    speed_periods = speed.reshape(-1,int(speed.size/4))  # reshape speed in baseline, low, mid, high injection time periods

    win_30 = fs*60*30 # 30 min window
    Lfp_B = lfp[0:win_30,:]  # base line period
    Lfp_L = lfp[win_30:2*win_30,:]  # low injection
    Lfp_M = lfp[2*win_30:3*win_30,:]  # mid injection
    Lfp_H = lfp[3*win_30:4*win_30,:]  # high injection

    speed_B = speed_periods[0,:]
    speed_L = speed_periods[1,:]
    speed_M = speed_periods[2,:]
    speed_H = speed_periods[3,:]

    print('min in each epoch: ',speed_B.size/fs/60)

    return Lfp_B, Lfp_L, Lfp_M, Lfp_H, speed_B, speed_L, speed_M, speed_H

# =============================================================================

''' Select only 1 min data at the time to speed up data processing, for both LFP and speed '''

def select_1min_data(Lfp_B, Lfp_L, Lfp_M, Lfp_H, speed_B, speed_L, speed_M, speed_H, current_min, offset, fs=2500):

    #### 3. Select 1 min of data for each Epoch, for both Lfp and speed
    ##### Select 1 every M channels for the Lfp
    # offset is the starting minute of epoch

    M = 1 # keep every M channel
    length = 1 # length of period to look at, i.e. 1 = 1 min

    start = (offset + current_min)*60*fs - 1    # start point in time points
    end = (offset + current_min + length)*60*fs - 1   # stop point in time points

    # Lfp trim: 1 minute, one every M channel and tranform memap to numpy array
    Lfp_B_min = np.array(Lfp_B[start:end,::M])  # base line period
    Lfp_L_min = np.array(Lfp_L[start:end,::M])  # low injection
    Lfp_M_min = np.array(Lfp_M[start:end,::M])  # mid injection
    Lfp_H_min = np.array(Lfp_H[start:end,::M])  # high injection

    # speed
    speed_B_min = speed_B[start:end]
    speed_L_min = speed_L[start:end]
    speed_M_min = speed_M[start:end]
    speed_H_min = speed_H[start:end]

    print('1 min data: Lfp shape {}, speed shape {}, length in sec: {}\n'.format(Lfp_B_min.shape, speed_B_min.shape, speed_B_min.size/fs))

    return Lfp_B_min, Lfp_L_min, Lfp_M_min, Lfp_H_min, speed_B_min, speed_L_min, speed_M_min, speed_H_min


# =============================================================================

''' 
Replace bad lfp channel with the nearest neighbor channel 
bad_id = bad channel id
next_id = nearest neighbor id
'''

def replace_bad_lfp_channel(Lfp_B_min, Lfp_L_min, Lfp_M_min, Lfp_H_min, bad_id, next_id):

    Lfp_B_min[:,bad_id] = Lfp_B_min[:,next_id]
    Lfp_L_min[:,bad_id] = Lfp_L_min[:,next_id]
    Lfp_M_min[:,bad_id] = Lfp_M_min[:,next_id]
    Lfp_H_min[:,bad_id] = Lfp_H_min[:,next_id]

    return Lfp_B_min, Lfp_L_min, Lfp_M_min, Lfp_H_min

# =============================================================================

"""
Average 4 channels in the Neuropixel together and keep only the average Lfp.
Channels averaged together are nearest neighbors in the x and y directions. 
They are in the same 2x2 block of channels.
The averaging is done for 1 min of data at the time to deal with excessive usage of RAM otherwise 
"""

def average_lfp_4_channels(Lfp_B_min,Lfp_L_min,Lfp_M_min,Lfp_H_min):


    print("Averaging lfp ...")

    if int(Lfp_B_min.shape[1]) % 2:
        sys.exit("Number of channels in the brain region considered is not EVEN! Check channel list!")

    # If the tot number of channels is not a multiple of 4, remove the extra channels
    if int(Lfp_B_min.shape[1] % 4):
        ch_extra = (Lfp_B_min.shape[1] % 4 ) # extra channels, which make the tot number not a multiple of 4
        Lfp_B_min = Lfp_B_min[:,:-ch_extra]
        Lfp_L_min = Lfp_L_min[:,:-ch_extra]
        Lfp_M_min = Lfp_M_min[:,:-ch_extra]
        Lfp_H_min = Lfp_H_min[:,:-ch_extra]
        print(f'-- {ch_extra} Lfp channels were removed in order to have the tot number of channels a multiple of 4, for the average in a 2x2 channel block\n')

    # baseline
    Lfp_B_min = Lfp_B_min - np.mean(Lfp_B_min, axis=1, keepdims=True) # reReferencing (CAR)
    Lfp_RS_B = Lfp_B_min.reshape(Lfp_B_min.shape[0],int(Lfp_B_min.shape[1]/4),4) # reshape Lfp, such that channels in the same 2x2 block are in the same colum dim
    Lfp_avg_B = Lfp_RS_B.mean(axis=2)
    # low injection
    Lfp_L_min = Lfp_L_min - np.mean(Lfp_L_min, axis=1, keepdims=True) # reReferencing (CAR)
    Lfp_RS_L = Lfp_L_min.reshape(Lfp_L_min.shape[0],int(Lfp_L_min.shape[1]/4),4) # reshape Lfp, such that channels in the same 2x2 block are in the same colum dim
    Lfp_avg_L = Lfp_RS_L.mean(axis=2)
    # mid injection
    Lfp_M_min = Lfp_M_min - np.mean(Lfp_M_min, axis=1, keepdims=True) # reReferencing (CAR)
    Lfp_RS_M = Lfp_M_min.reshape(Lfp_M_min.shape[0],int(Lfp_M_min.shape[1]/4),4) # reshape Lfp, such that channels in the same 2x2 block are in the same colum dim
    Lfp_avg_M = Lfp_RS_M.mean(axis=2)
    # high injection
    Lfp_H_min = Lfp_H_min - np.mean(Lfp_H_min, axis=1, keepdims=True) # reReferencing (CAR)
    Lfp_RS_H = Lfp_H_min.reshape(Lfp_H_min.shape[0],int(Lfp_H_min.shape[1]/4),4) # reshape Lfp, such that channels in the same 2x2 block are in the same colum dim
    Lfp_avg_H = Lfp_RS_H.mean(axis=2)

    return Lfp_avg_B, Lfp_avg_L, Lfp_avg_M, Lfp_avg_H

# =============================================================================

"""
Average 2 channels in the Neuropixel array which are on the same depth, same y 
i.e. average consecutive channels in the Lfp map, for each epoch separately, for 1 min of data at the time 
"""
def average_lfp_same_depth(Lfp_B_min,Lfp_L_min,Lfp_M_min,Lfp_H_min):

    print("Averaging lfp ...")

    if int(Lfp_B_min.shape[1]) % 2:
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
    Lfp_RS_H = Lfp_H_min.reshape(Lfp_H_min.shape[0],int(Lfp_H_min.shape[1]/2),2) # reshape Lfp, such that channels on the same depth are into adjacent columns
    Lfp_avg_H = Lfp_RS_H.mean(axis=2)

    return Lfp_avg_B, Lfp_avg_L, Lfp_avg_M, Lfp_avg_H


# =============================================================================

"""" Subsample LFPs and speed to 1250 Hz """

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

""" Stack all the 20 min of recordings for each epoch together, (n_min, T, n_ch) """

def stack_lfp_1min_all_trials(lfp_B_epoch,lfp_L_epoch,lfp_M_epoch,lfp_H_epoch, lfp_dec_B, lfp_dec_L, lfp_dec_M, lfp_dec_H):


    lfp_B_epoch.append(lfp_dec_B) # add one minute lfp: trial num x trial length, for each channel
    lfp_L_epoch.append(lfp_dec_L)
    lfp_M_epoch.append(lfp_dec_M)
    lfp_H_epoch.append(lfp_dec_H)



    return lfp_B_epoch, lfp_L_epoch, lfp_M_epoch, lfp_H_epoch

# =============================================================================


""" 
Create a mask for the LFP trials with artifacts 

input: 
lfp_dec: lfp (subsampled)
win: time length of the window we are looking at in time points, i.e. 1 sec = 1250 points
std_th: std threshold
period: string for baseline, low, etc..

"""
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

"""
Create a speed mask. 

True if speed < threshold at any point in a given time window of length win
False if there is at least one value above threshold in the given time window
Input: 
speed_up: upsampled speed
win: amplitude in data points of the win we are looking at for thresholding
thresh: speed threshold for speed
period: string for baseline, low, mid, high injection 

"""

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

"""
Create mask for:
    1. Lfp artifacts
    2. speed (low and high speed)
    3. Combine them together

win = window length, e.g. 1250 = 1 sec
th = speed threshold 

Input: 1 min Lfp in the form: time x channel
Output: mask for low (high) speed trial together with LFP artifacts trials, for each epoch

"""

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

""" Stack mask for low and high speed relative to a given minute, in order to have masks for the whole 20 min period """

def stack_mask_1min(mask_B_low, mask_L_low, mask_M_low, mask_H_low,
                     mask_B_high, mask_L_high, mask_M_high, mask_H_high,
                     tot_mask_B_low_s, tot_mask_L_low_s, tot_mask_M_low_s, tot_mask_H_low_s,
                     tot_mask_B_high_s, tot_mask_L_high_s, tot_mask_M_high_s, tot_mask_H_high_s):


    mask_B_low.append(tot_mask_B_low_s) # add one minute lfp: trial num x trial length, for each channel
    mask_L_low.append(tot_mask_L_low_s)
    mask_M_low.append(tot_mask_M_low_s)
    mask_H_low.append(tot_mask_H_low_s)

    mask_B_high.append(tot_mask_B_high_s) # add one minute lfp: trial num x trial length, for each channel
    mask_L_high.append(tot_mask_L_high_s)
    mask_M_high.append(tot_mask_M_high_s)
    mask_H_high.append(tot_mask_H_high_s)



    return mask_B_low, mask_L_low, mask_M_low, mask_H_low, mask_B_high, mask_L_high, mask_M_high, mask_H_high

# =============================================================================


"""
Keep only LFP good trials: those without artifacts
Method: mask LFP with the mask-good-trials
"""

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

""" 
For each channel, stack all the min recordings together. For each epoch.
Output: (min, T, channel) 
"""
def stack_lfp_1min(lfp_B_epoch,lfp_L_epoch,lfp_M_epoch,lfp_H_epoch,lfp_B_list,lfp_L_list,lfp_M_list,lfp_H_list):

    nch = len(lfp_B_list)
    for ch in range(nch):
        lfp_B_epoch[ch].append(lfp_B_list[ch]) # add one minute lfp: trial num x trial length, for each channel
        lfp_L_epoch[ch].append(lfp_L_list[ch])
        lfp_M_epoch[ch].append(lfp_M_list[ch])
        lfp_H_epoch[ch].append(lfp_H_list[ch])


    return lfp_B_epoch, lfp_L_epoch, lfp_M_epoch, lfp_H_epoch




# =============================================================================

# Save low and high speed LFP trials (only good trials, i.e. without artifact),
# shaped into a 4D array nch, min id, id trial, length trial

def save_matlab_files(rec, sess, brain_reg,
                      lfp_B_ep_low_s, lfp_L_ep_low_s, lfp_M_ep_low_s,lfp_H_ep_low_s,
                      lfp_B_ep_high_s, lfp_L_ep_high_s, lfp_M_ep_high_s, lfp_H_ep_high_s, save_var):

    main_dir = r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC'
    path = rec[sess][brain_reg]

    dir_sess = path.split('\\')[-3] # path for session directory
    full_dir_path = os.path.join(main_dir,dir_sess,'LFPs_and_masks')

    if not os.path.exists(full_dir_path):
        os.makedirs(full_dir_path)

    # =============================================================================
    #     Low Speed Lfp
    # =============================================================================

    # file path to save in matlab
    out_file = os.path.join(full_dir_path, "lfp_epoch_low_speed_{}.mat".format(save_var))

    mat_lfp = {'B': lfp_B_ep_low_s,
               'L': lfp_L_ep_low_s,
               'M': lfp_M_ep_low_s,
               'H': lfp_H_ep_low_s}

    data_lfp = {'lfp': mat_lfp}

    # save lfp for each epoch in matlab format
    savemat(out_file, data_lfp)


    # =============================================================================
    #     High Speed Lfp
    # =========================================================================

    # file path to save in matlab
    out_file = os.path.join(full_dir_path, "lfp_epoch_high_speed_{}.mat".format(save_var))

    mat_lfp = {'B': lfp_B_ep_high_s,
               'L': lfp_L_ep_high_s,
               'M': lfp_M_ep_high_s,
               'H': lfp_H_ep_high_s}

    data_lfp = {'lfp': mat_lfp}

    # save lfp for each epoch in matlab format
    savemat(out_file, data_lfp)



# =============================================================================

# Save the whole Lfp in the form: min x time x channel
# Save masks for low speed-no artifact and high speed-no artifcats

# def save_matlab_files_all_lfps(rec,sess,brain_reg, lfp_B_ep, lfp_L_ep, lfp_M_ep, lfp_H_ep,
#                                tot_mask_B_low_s, tot_mask_L_low_s, tot_mask_M_low_s, tot_mask_H_low_s,
#                                tot_mask_B_high_s, tot_mask_L_high_s, tot_mask_M_high_s, tot_mask_H_high_s):

#     main_dir = r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC'
#     path = rec[sess][brain_reg]

#     dir_sess = path.split('\\')[-3] # path for session directory
#     full_dir_path = os.path.join(main_dir,dir_sess)

#     if not os.path.exists(full_dir_path):
#         os.makedirs(full_dir_path)


#     # =============================================================================
#     #     Lfp all trials
#     # =============================================================================

#     # file path to save in matlab
#     out_file = os.path.join(full_dir_path, "lfp_epoch_high_speed_CSD.mat")

#     mat_lfp = {'B': lfp_B_ep_high_s,
#                'L': lfp_L_ep_high_s,
#                'M': lfp_M_ep_high_s,
#                'H': lfp_H_ep_high_s}

#     data_lfp = {'lfp': mat_lfp}

#     # save lfp for each epoch in matlab format
#     savemat(out_file, data_lfp)




# =============================================================================

# Save the whole Lfp in the form: min x time x channel
# Save masks for low speed-no artifact and high speed-no artifcats

def save_matlab_files_all_lfps(rec,sess,brain_reg, lfp_B_ep, lfp_L_ep, lfp_M_ep, lfp_H_ep,
                               tot_mask_B_low_s, tot_mask_L_low_s, tot_mask_M_low_s, tot_mask_H_low_s,
                               tot_mask_B_high_s, tot_mask_L_high_s, tot_mask_M_high_s, tot_mask_H_high_s, save_var):

    main_dir = r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC'
    path = rec[sess][brain_reg]

    dir_sess = path.split('\\')[-3] # path for session directory
    full_dir_path = os.path.join(main_dir,dir_sess,'LFPs_and_masks')

    if not os.path.exists(full_dir_path):
        os.makedirs(full_dir_path)


    # =============================================================================
    #     Lfp all trials
    # =============================================================================

    # file path to save in matlab
    out_file = os.path.join(full_dir_path, "lfp_epoch_all_trials_{}.mat".format(save_var))

    mat_lfp = {'B': np.array(lfp_B_ep),
               'L': np.array(lfp_L_ep),
               'M': np.array(lfp_M_ep),
               'H': np.array(lfp_H_ep)}

    data_lfp_all = {'lfp_all': mat_lfp}

    # save lfp for each epoch in matlab format
    savemat(out_file, data_lfp_all)


    # =============================================================================
    #     Mask low/high speed - no artifact
    # =============================================================================

    # file path to save Masks, low speed-no artifact in matlab

    out_file = os.path.join(full_dir_path, "mask_low_high_speed_{}.mat".format(save_var))


    # create dictionaries to save in matlab
    mat_mask = {'B_low': tot_mask_B_low_s,
                'L_low': tot_mask_L_low_s,
                'M_low': tot_mask_M_low_s,
                'H_low': tot_mask_H_low_s,
                'B_high': tot_mask_B_high_s,
                'L_high': tot_mask_L_high_s,
                'M_high': tot_mask_M_high_s,
                'H_high': tot_mask_H_high_s}


    # save lfp for each epoch in matlab format
    data_mask = {'mask': mat_mask}

    # save lfp for each epoch in matlab format
    savemat(out_file, data_mask)



# =============================================================================



"""
FOR HIGH LARGE MATRICES: i.e. All lfp channels 

Save low and high speed LFP trials (only good trials, i.e. without artifact),
shaped into a 4D array nch, min id, id trial, length trial
"""


def save_h5_files(rec, sess, brain_reg, lfp_B_ep_low_s, lfp_L_ep_low_s, lfp_M_ep_low_s, lfp_H_ep_low_s,
                  lfp_B_ep_high_s, lfp_L_ep_high_s, lfp_M_ep_high_s, lfp_H_ep_high_s):

    main_dir = r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC'
    path = rec[sess][brain_reg]

    dir_sess = path.split('\\')[-3] # path for session directory
    full_dir_path = os.path.join(main_dir, dir_sess)

    if not os.path.exists(full_dir_path):
        os.makedirs(full_dir_path)

    # =============================================================================
    #     Low Speed Lfp
    # =============================================================================

    # file path to save in hdf5
    lfp_low_speed_file = os.path.join(full_dir_path, "lfp_epoch_low_speed_all_ch.h5")

    with h5py.File(lfp_low_speed_file, 'w') as f:
        f.create_dataset('B', data=np.array(lfp_B_ep_low_s, dtype=np.float64))
        f.create_dataset('L', data=np.array(lfp_L_ep_low_s, dtype=np.float64))
        f.create_dataset('M', data=np.array(lfp_M_ep_low_s, dtype=np.float64))
        f.create_dataset('H', data=np.array(lfp_H_ep_low_s, dtype=np.float64))

    # =============================================================================
    #     High Speed Lfp
    # =============================================================================

    # file path to save in hdf5
    lfp_high_speed_file = os.path.join(full_dir_path, "lfp_epoch_high_speed_all_ch.h5")

    with h5py.File(lfp_high_speed_file, 'w') as f:
        f.create_dataset('B', data=np.array(lfp_B_ep_high_s, dtype=np.float64))
        f.create_dataset('L', data=np.array(lfp_L_ep_high_s, dtype=np.float64))
        f.create_dataset('M', data=np.array(lfp_M_ep_high_s, dtype=np.float64))
        f.create_dataset('H', data=np.array(lfp_H_ep_high_s, dtype=np.float64))

# Usage example:
# save_h5_files(rec, sess, brain_reg, lfp_B_ep_low_s, lfp_L_ep_low_s, lfp_M_ep_low_s, lfp_H_ep_low_s,
#               lfp_B_ep_high_s, lfp_L_ep_high_s, lfp_M_ep_high_s, lfp_H_ep_high_s)


# =============================================================================

def save_hdf5_files(rec, sess, brain_reg, lfp_B_ep_low_s, lfp_L_ep_low_s, lfp_M_ep_low_s, lfp_H_ep_low_s,
                    lfp_B_ep_high_s, lfp_L_ep_high_s, lfp_M_ep_high_s, lfp_H_ep_high_s):

    main_dir = r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC'
    path = rec[sess][brain_reg]

    dir_sess = path.split('\\')[-3] # path for session directory
    full_dir_path = os.path.join(main_dir, dir_sess)

    if not os.path.exists(full_dir_path):
        os.makedirs(full_dir_path)

    def save_list_of_arrays(h5file, group_name, list_of_arrays):
        group = h5file.create_group(group_name)
        for i, arr in enumerate(list_of_arrays):
            group.create_dataset(str(i), data=arr)

    # =============================================================================
    #     Low Speed Lfp
    # =============================================================================

    lfp_low_speed_file = os.path.join(full_dir_path, "lfp_epoch_low_speed_all_ch.h5")
    with h5py.File(lfp_low_speed_file, 'w') as h5f:
        save_list_of_arrays(h5f, 'B', lfp_B_ep_low_s)
        save_list_of_arrays(h5f, 'L', lfp_L_ep_low_s)
        save_list_of_arrays(h5f, 'M', lfp_M_ep_low_s)
        save_list_of_arrays(h5f, 'H', lfp_H_ep_low_s)

    # =============================================================================
    #     High Speed Lfp
    # =============================================================================

    lfp_high_speed_file = os.path.join(full_dir_path, "lfp_epoch_high_speed_all_ch.h5")
    with h5py.File(lfp_high_speed_file, 'w') as h5f:
        save_list_of_arrays(h5f, 'B', lfp_B_ep_high_s)
        save_list_of_arrays(h5f, 'L', lfp_L_ep_high_s)
        save_list_of_arrays(h5f, 'M', lfp_M_ep_high_s)
        save_list_of_arrays(h5f, 'H', lfp_H_ep_high_s)


# ============================================================================

"""
FOR HIGH LARGE MATRICES: i.e. All lfp channels 

Save files lfps (all trials) and masks to be opened in matlab. 
This function uses h5py in order to save large files
"""

def save_matlab_files_all_lfps_h5(rec, sess, brain_reg, lfp_B_ep, lfp_L_ep, lfp_M_ep, lfp_H_ep,
                               tot_mask_B_low_s, tot_mask_L_low_s, tot_mask_M_low_s, tot_mask_H_low_s,
                               tot_mask_B_high_s, tot_mask_L_high_s, tot_mask_M_high_s, tot_mask_H_high_s):

    main_dir = r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC'
    path = rec[sess][brain_reg]

    dir_sess = path.split('\\')[-3]  # path for session directory
    full_dir_path = os.path.join(main_dir, dir_sess)

    if not os.path.exists(full_dir_path):
        os.makedirs(full_dir_path)

    # Save LFPs
    lfp_file = os.path.join(full_dir_path, "lfp_epoch_all_trials_all_ch.h5")
    with h5py.File(lfp_file, 'w') as f:
        f.create_dataset('B', data=np.array(lfp_B_ep, dtype=np.float64))
        f.create_dataset('L', data=np.array(lfp_L_ep, dtype=np.float64))
        f.create_dataset('M', data=np.array(lfp_M_ep, dtype=np.float64))
        f.create_dataset('H', data=np.array(lfp_H_ep, dtype=np.float64))

    # Save Masks
    mask_file = os.path.join(full_dir_path, "mask_low_high_speed_all_ch.h5")
    with h5py.File(mask_file, 'w') as f:
        f.create_dataset('B_low', data=np.array(tot_mask_B_low_s, dtype=np.bool_))
        f.create_dataset('L_low', data=np.array(tot_mask_L_low_s, dtype=np.bool_))
        f.create_dataset('M_low', data=np.array(tot_mask_M_low_s, dtype=np.bool_))
        f.create_dataset('H_low', data=np.array(tot_mask_H_low_s, dtype=np.bool_))
        f.create_dataset('B_high', data=np.array(tot_mask_B_high_s, dtype=np.bool_))
        f.create_dataset('L_high', data=np.array(tot_mask_L_high_s, dtype=np.bool_))
        f.create_dataset('M_high', data=np.array(tot_mask_M_high_s, dtype=np.bool_))
        f.create_dataset('H_high', data=np.array(tot_mask_H_high_s, dtype=np.bool_))

# =============================================================================


"""
Convert a list of lists of arrays (channels, minutes, M, T) into a 3D array (nch, N, T).

:param lfp_data: A list (channels) of lists (minutes) of arrays (M x T) with M variable.
:return: A 3D NumPy array with shape (nch, N, T).
"""


def stack_lfp_data_to_array(lfp_data):

    channel_arrays = []

    # Calculate the total size N for each channel
    for channel_data in lfp_data:
        # Concatenate all arrays from all minutes for the current channel
        concatenated_data = np.concatenate([array for minute_data in channel_data for array in minute_data], axis=0)
        channel_arrays.append(concatenated_data)

    # Stack along a new axis to get a 3D array (nch, N, T)
    return np.stack(channel_arrays, axis=0)

# Example usage:
# Assuming `lfp_data` is your list of lists of arrays with the structure (ch, min, M, T)
# lfp_array = stack_lfp_data_to_array(lfp_data)
# Now `lfp_array` is a 3D NumPy array with shape (nch, N, T)



# =============================================================================


""" -------------------------------------------------------
Save list of electrodes in CA1
"""

def save_list_channel_CA1(idx_cell_hpc, rec, sess, brain_reg):

    main_dir = r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC'
    path = rec[sess][brain_reg]

    dir_sess = path.split('\\')[-3]     # path for session directory
    full_dir_path = os.path.join(main_dir, dir_sess, 'freq_phase_matrices')

    np.save(os.path.join(full_dir_path,'idx_cell_hpc.npy'),idx_cell_hpc)


def load_list_channel_CA1(rec, sess, brain_reg):

    main_dir = r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC'
    path = rec[sess][brain_reg]

    dir_sess = path.split('\\')[-3]     # path for session directory
    full_dir_path = os.path.join(main_dir, dir_sess, 'freq_phase_matrices')

    idx_cell_hpc = np.load(os.path.join(full_dir_path,'idx_cell_hpc.npy'))

    return idx_cell_hpc
