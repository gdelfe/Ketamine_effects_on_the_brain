"""
Spyder Editor

@ Gino Del Ferraro, Fenton lab, Oct 2023
"""

import numpy as np

from scipy.signal import butter, lfilter, freqz, iirnotch, welch, filtfilt, iirfilter, iirfilter 
from scipy.interpolate import interp1d
from scipy.stats import zscore
import sys

# =============================================================================
# LFP filtering functions 
# =============================================================================

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


def bandstop_filter(data, f_low, f_high, fs=2500, order=4):
    nyq = 0.5 * fs
    low = f_low / nyq
    high = f_high / nyq

    b, a = iirfilter(N=order, Wn=[low, high], btype='bandstop', ftype='butter')
    
    filtered_data = lfilter(b, a, data)
    
    return filtered_data


# =============================================================================
# Speed-related functions 
# =============================================================================


def upsample_speed(speed, LFP, sess, LFP_rate = 2500, speed_rate = 100):
    
    # time length for speed and LFP variables
    speed_T = np.linspace(0, len(speed[sess])/ speed_rate, len(speed[sess]))
    LFP_T = np.linspace(0, len(LFP[:,0])/ LFP_rate, len(LFP[:,0])) # all channels have same length

    # interpolate speed variable based on LFP
    interpolator = interp1d(speed_T, speed[sess], kind = 'linear', fill_value="extrapolate")
    speed_upsampled = interpolator(LFP_T) 
    
    # check if upsample done correctly
    if (speed_upsampled.size - LFP.shape[0]) !=0:
        sys.exit("Speed upsampled size and Lfp size are not the same!")
    print('speed upsampled shape {}, Lfp shape {}'.format(speed_upsampled.shape, LFP.shape))
    
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


# =============================================================================

# scale Lfp into mV and band pass it at [1,300] Hz by using a non-causal filter
# Apply notch filter at 60 Hz to remove power line distortion

def filter_lfp_in_each_epoch(Lfp_B_min,Lfp_L_min,Lfp_M_min,Lfp_H_min,gain,qband):

    print('Filtering Lfp ...')
    
    # baseline
    lfp_scaled_B = Lfp_B_min*gain*1e6 # scale in uV
    lfp_filt_B_bp = bandpass_filter(lfp_scaled_B, lowcut = 1, highcut = 300, fs=2500, order=5) # band pass filter at 1 Hz and 300 Hz
    lfp_filt_B = notch_filter(lfp_filt_B_bp, 60, fs=2500, Q=qband)
    
    # low injection 
    lfp_scaled_L = Lfp_L_min*gain*1e6 # scale in uV
    lfp_filt_L_bp = bandpass_filter(lfp_scaled_L, lowcut = 1, highcut = 300, fs=2500, order=5) # band pass filter at 1 Hz and 300 Hz
    lfp_filt_L = notch_filter(lfp_filt_L_bp, 60, fs=2500, Q=qband)
    
    # mid injection 
    lfp_scaled_M = Lfp_M_min*gain*1e6 # scale in uV
    lfp_filt_M_bp = bandpass_filter(lfp_scaled_M, lowcut = 1, highcut = 300, fs=2500, order=5) # band pass filter at 1 Hz and 300 Hz
    lfp_filt_M = notch_filter(lfp_filt_M_bp, 60, fs=2500, Q=qband)
    
    # high injection 
    lfp_scaled_H = Lfp_H_min*gain*1e6 # scale in uV
    lfp_filt_H_bp = bandpass_filter(lfp_scaled_H, lowcut = 1, highcut = 300, fs=2500, order=5) # band pass filter at 1 Hz and 300 Hz
    lfp_filt_H = notch_filter(lfp_filt_H_bp, 60, fs=2500, Q=qband)

    return lfp_filt_B, lfp_filt_L, lfp_filt_M, lfp_filt_H

# =============================================================================

# Compute Current Source Density.
# Output: CSD and CSD filtered 

def compute_iCSD(Lfp):
    
    
    '''iCSD toolbox demonstration script'''
    import matplotlib.pyplot as plt
    import numpy as np
    import icsd
    from scipy import io
    import neo
    import quantities as pq
    
    lastdefinition = pq.A / pq.V  # This line initializes lastdefinition to the base unit for Siemens.
    
    #patch quantities with the SI unit Siemens if it does not exist
    for symbol, prefix, definition, u_symbol in zip(
        ['siemens', 'S', 'mS', 'uS', 'nS', 'pS'],
        ['', '', 'milli', 'micro', 'nano', 'pico'],
        [pq.A/pq.V, pq.A/pq.V, 'S', 'mS', 'uS', 'nS'],
        [None, None, None, None, u'µS', None]):
        if type(definition) is str:
            definition = lastdefinition / 1000
        if not hasattr(pq, symbol):
            setattr(pq, symbol, pq.UnitQuantity(
                prefix + 'siemens',
                definition,
                symbol=symbol,
                u_symbol=u_symbol))
        lastdefinition = definition
    
    
    #prepare lfp data for use, by changing the units to SI and append quantities,
    #along with electrode geometry, conductivities and assumed source geometry
    lfp_data = Lfp.T * 1E-6 * pq.V        # [uV] -> [V]
    z_data = np.linspace(20E-6,20E-6*Lfp.shape[1],Lfp.shape[1]) * pq.m  # [m]
    diam = 500E-6 * pq.m                              # [m]
    h = 20E-6 * pq.m                                 # [m]
    sigma = 0.3 * pq.S / pq.m                         # [S/m] or [1/(ohm*m)]
    sigma_top = 0.3 * pq.S / pq.m                     # [S/m] or [1/(ohm*m)]
    
    
    spline_input = {
    'lfp' : lfp_data,
    'coord_electrode' : z_data,
    'diam' : diam,
    'sigma' : sigma,
    'sigma_top' : sigma,
    'num_steps' : Lfp.shape[1],      # Spatial CSD upsampling to N steps
    'tol' : 1E-12,
    'f_type' : 'gaussian',
    'f_order' : (20, 5),
    }
    
    #Create the different CSD-method class instances. We use the class methods
    #get_csd() and filter_csd() below to get the raw and spatially filtered
    #versions of the current-source density estimates.
    csd_dict = dict(spline_icsd = icsd.SplineiCSD(**spline_input))
    
    #plot
    for method, csd_obj in list(csd_dict.items()):
        fig, axes = plt.subplots(3,1, figsize=(8,8))
        
        #plot LFP signal
        ax = axes[0]
        im = ax.imshow(np.array(lfp_data), origin='upper', vmin=-abs(lfp_data).max(), \
                  vmax=abs(lfp_data).max(), cmap='jet_r', interpolation='nearest')
        ax.axis(ax.axis('tight'))
        cb = plt.colorbar(im, ax=ax)
        cb.set_label('LFP (%s)' % lfp_data.dimensionality.string)
        ax.set_xticklabels([])
        ax.set_title('LFP')
        ax.set_ylabel('ch #')
        
        #plot raw csd estimate
        csd = csd_obj.get_csd()
        ax = axes[1]
        im = ax.imshow(np.array(csd), origin='upper', vmin=-abs(csd).max(), \
              vmax=abs(csd).max(), cmap='jet_r', interpolation='nearest')
        ax.axis(ax.axis('tight'))
        ax.set_title(csd_obj.name)
        cb = plt.colorbar(im, ax=ax)
        cb.set_label('CSD (%s)' % csd.dimensionality.string)
        ax.set_xticklabels([])
        ax.set_ylabel('ch #')
        
        #plot spatially filtered csd estimate
        ax = axes[2]
        csd_fil = csd_obj.filter_csd(csd)
        im = ax.imshow(np.array(csd_fil), origin='upper', vmin=-abs(csd_fil).max(), \
              vmax=abs(csd_fil).max(), cmap='jet_r', interpolation='nearest')
        ax.axis(ax.axis('tight'))
        ax.set_title(csd_obj.name + ', filtered')
        cb = plt.colorbar(im, ax=ax)
        cb.set_label('CSD (%s)' % csd_fil.dimensionality.string)
        ax.set_ylabel('ch #')
        ax.set_xlabel('timestep')
        
    return csd.T, csd_fil.T 



