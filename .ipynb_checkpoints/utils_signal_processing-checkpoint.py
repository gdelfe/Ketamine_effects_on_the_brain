import numpy as np

from scipy.signal import butter, lfilter, freqz, iirnotch, welch, filtfilt 
from scipy.interpolate import interp1d
from scipy.stats import zscore
import sys

# -----------------------------------------------------
# LFP filtering functions 
# -----------------------------------------------------

def bandpass_filter(data, lowcut, highcut, fs=2500, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')

    filtered_data = np.empty_like(data)

    # Loop through each channel and apply the filter
    for i in range(data.shape[1]):
        filtered_data[:, i] = filtfilt(b, a, data[:, i])

    return filtered_data

# -----------------------------------------------------

def highpass_filter(data, highcut, fs=2500, order=5):
      
    nyq = 0.5 * fs
    
    # high-pass filter
    high = highcut / nyq
    b, a = butter(order, high, btype='high')
    
    filtered_data = np.empty_like(data)

    # Loop through each channel and apply the filter
    for i in range(data.shape[1]):
        filtered_data[:, i] = filtfilt(b, a, data[:, i])
        
    return filtered_data

# -----------------------------------------------------

def notch_filter(data, notch_freq, fs=2500, Q=30):
    b, a = iirnotch(notch_freq, Q, fs)
    filtered_data = np.empty_like(data)

    # Loop through each channel and apply the filter
    for i in range(data.shape[1]):
        filtered_data[:, i] = lfilter(b, a, data[:, i])

    return filtered_data


# -----------------------------------------------------
# Speed-related functions 
# -----------------------------------------------------

def upsample_speed(speed, LFP, sess, LFP_rate = 2500, speed_rate = 100):
    
    # time length for speed and LFP variables
    speed_T = np.linspace(0, len(speed[sess])/ speed_rate, len(speed[sess]))
    LFP_T = np.linspace(0, len(LFP[:,0])/ LFP_rate, len(LFP[:,0])) # all channels have same length

    # interpolate speed variable based on LFP
    interpolator = interp1d(speed_T, speed[sess], kind = 'linear', fill_value="extrapolate")
    speed_upsampled = interpolator(LFP_T) 
    
    return speed_upsampled

# -----------------------------------------------------
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


# -----------------------------------------------------
# Lfp artifact identification and masking 
# -----------------------------------------------------

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