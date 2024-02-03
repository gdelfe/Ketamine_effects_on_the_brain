# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 15:01:42 2023

@author: Gino Del Ferraro

"""

import numpy as np
import os
import pickle
from utilities_ketamine_analysis_v8 import *
from utils_signal_processing import *
from utils_plotting import *
from utils_general import *
from utils_phase_frequency import *
from scipy.signal import morlet, welch, convolve
from scipy.io import loadmat 
import matplotlib.gridspec as gridspec
import time 

sess = 2
hpc_start = 75 
hpc_end = 131 

name_lfp_file = "lfp_epoch_all_trials_CSD_x_avg.mat"
mask_low_and_high_speed_file = "mask_low_high_speed_CSD_x_avg.mat"

# --------------

print('N of channels in CA1', hpc_end - hpc_start)


""" load Path for ephys """


out_path = Path(r'C:\Users\fentonlab\Desktop\Gino\HPC_channel_info\\')

sync_PATH = os.path.join(out_path,'sync_path_matched.pkl')
nt_PATH = os.path.join(out_path,'nt_path_matched.pkl')
ksort_PATH = os.path.join(out_path,'ksort_HPC_path_matched.pkl')
ts_ephys_PATH = os.path.join(out_path,'ts_ephys_aln.pkl')

with open(sync_PATH,'rb') as file:
    sync_path_matched = pickle.load(file)
    
with open(nt_PATH,'rb') as file:
    nt_path_matched= pickle.load(file)
    
with open(ksort_PATH,'rb') as file:
    ksort_HPC_path_matched = pickle.load(file)
    
with open(ts_ephys_PATH,'rb') as file:
    ts_ephys_aln = pickle.load(file)
    
#%%    
    
""" Load spike files"""

start_time = time.time()
# load HPC file names and store them in rec spikes at 250 Hz
path_spike = Path(r'\\sshfs.r\gino@monk.cns.nyu.edu\f\fentonlab\RAWDATA\NeuroPix\spk_ketamine\HPC\\')
rec_name = 'spk_250Hz_' + str(sess) + '.file'
spk_file_name = os.path.join(path_spike , rec_name)

with open(spk_file_name,'rb') as f:
    spikes_binned_all_HPC = pickle.load(f)
    
time_taken = time.time() - start_time
print(f"Time of execution {time_taken} in sec")


#%%
# check dimension of spikes and behaviour (at 100 Hz)
spikes_binned_all_HPC[0].shape[1]

spk = spikes_binned_all_HPC
 

"""Load LFP for a given session -- Only CA1"""
    
binFullPath = r'C:\Users\fentonlab\Desktop\Gino\LFPs'
main_dir = r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC'
HPC_path_file = os.path.join(r'C:\Users\fentonlab\Desktop\Gino\LFPs','HPC_lfp_paths.file')
PFC_path_file = os.path.join(r'C:\Users\fentonlab\Desktop\Gino\LFPs','PFC_lfp_paths.file')

# ====== Load Lfp and speed data for a specific recording and brain area 
rec = load_rec_path(binFullPath,HPC_path_file,PFC_path_file,"HPC",sess)
path = rec[sess]['HPC']
dir_sess = path.split('\\')[-3] # path for session directory
full_dir_path = os.path.join(main_dir,dir_sess)
out_file = os.path.join(full_dir_path, name_lfp_file)
lfp = loadmat(out_file)


# assign Lfp varibales
lfp_B = lfp['lfp_all']['B'][0][0]
lfp_L = lfp['lfp_all']['L'][0][0]
lfp_M = lfp['lfp_all']['M'][0][0]
lfp_H = lfp['lfp_all']['H'][0][0]

lfp_B.shape

#### Downsample Lfp at 250 Hz (input at 1250 Hz)
LfpSB = lfp_B[:,::5,:] # original Lfp is sampled at 1250 Hz 
LfpSL = lfp_L[:,::5,:] 
LfpSM = lfp_M[:,::5,:] 
LfpSH = lfp_H[:,::5,:] 

LfpSB.shape


### Load low/high speed masks for a given session
in_file = os.path.join(full_dir_path, mask_low_and_high_speed_file)
mask = loadmat(in_file)

mask_B = mask['mask']['B_low'][0][0]
mask_L = mask['mask']['L_low'][0][0]
mask_M = mask['mask']['M_low'][0][0]
mask_H = mask['mask']['H_low'][0][0]


### spike-cell index for cells in HPC

# spike time of cell 0
spk[0][0][0:1000]
# mapping array between cell number and LFP channel
spk[2][:]


condition = (spk[2] >= hpc_start*2) & (spk[2] <= hpc_end*2) # get map of cells in CA1
idx_cell_hpc = np.where(condition)[0] # get index of cells in CA1
print('lfp cell index', idx_cell_hpc)



condition = (spk[2] >= hpc_start*2) & (spk[2] <= hpc_end*2) # get map of cells in CA1
idx_cell_hpc = np.where(condition)[0] # get index of cells in CA1
print("Cell index in full labeling:  \n", idx_cell_hpc) # print cell index 
print("Lfp channel associated to cell in full labeling:  \n", spk[2][idx_cell_hpc]) # print cell LFP channel 

#### Frequencies in logscale for spike-LFP analysis
freq = np.logspace(0.32, 2.48, num=70, endpoint=True)
freq


#%%

"""
**** CHOOSE YOUR CHANNEL LIST HERE ****
Channel indexes of the different CA1 strata (oriens, pyramidale, radiatum, loc mol, 
"""
lfp_ch_idx = list(range(6,11)) # CA1 layer 
# baseline 
lfp_B_CA1 = lfp_B[:,:,lfp_ch_idx]
mask_B_CA1 = mask_B[:,:,lfp_ch_idx]
# low dosage
lfp_L_CA1 = lfp_L[:,:,lfp_ch_idx]
mask_L_CA1 = mask_L[:,:,lfp_ch_idx]
# mid dosage
lfp_M_CA1 = lfp_M[:,:,lfp_ch_idx]
mask_M_CA1 = mask_M[:,:,lfp_ch_idx]
# high dosage
lfp_H_CA1 = lfp_H[:,:,lfp_ch_idx]
mask_H_CA1 = mask_H[:,:,lfp_ch_idx]

mask_B_CA1.shape, lfp_B_CA1.shape


"""
PHASE EXTRACTION FROM VARIOUS LFP CHANNELS in HPC layer 
"""
fs = 1250
n_cycles = 10

start_time = time.time()
# mask for low speed LFP, same sampling rate than the LFP mask_B_CA1,fs,axis=1) # 0 / 1 mask, same length of LFP --- used for calculations, only CA1
rep_mask_B_CA1 = np.repeat(mask_B_CA1,fs,axis=1)
rep_mask_L_CA1 = np.repeat(mask_L_CA1,fs,axis=1) 
rep_mask_M_CA1 = np.repeat(mask_M_CA1,fs,axis=1) 
rep_mask_H_CA1 = np.repeat(mask_H_CA1,fs,axis=1) 

# Phase and Convolution (envelope) computation for low speed trial only

phase_B_up, conv_B_up = phase_extraction_lfp_low_speed(lfp_B_CA1, rep_mask_B_CA1, freq, n_cycles, fs)
phase_L_up, conv_L_up = phase_extraction_lfp_low_speed(lfp_L_CA1, rep_mask_L_CA1, freq, n_cycles, fs)
phase_M_up, conv_M_up = phase_extraction_lfp_low_speed(lfp_M_CA1, rep_mask_M_CA1, freq, n_cycles, fs)
phase_H_up, conv_H_up = phase_extraction_lfp_low_speed(lfp_H_CA1, rep_mask_H_CA1, freq, n_cycles, fs)
time_taken = time.time() - start_time

print(f"Time of execution {time_taken} in sec")
print(f"Time of execution {time_taken/60} in min")

#%%
### 4.C Downsample phase to 250 Hz -- original phase at 1250 Hz

# baseline
phase_B = phase_B_up[:,::5,:,:]
mask_Bd = rep_mask_B_CA1[:,::5,:] 
# low dose
phase_L = phase_L_up[:,::5,:,:]
mask_Ld = rep_mask_L_CA1[:,::5,:] 
# mid dose
phase_M = phase_M_up[:,::5,:,:]
mask_Md = rep_mask_M_CA1[:,::5,:] 
# high dose 
phase_H = phase_H_up[:,::5,:,:]
mask_Hd = rep_mask_H_CA1[:,::5,:] 
phase_B.shape

## D. Phase-Frequency plots for each epoch
### 1.D Reshape phase in (entire period length, ch, frequency), period length is 20 min!

phase_BR = phase_B.reshape(-1,phase_B.shape[2],phase_B.shape[3])   # (T, ch, freq)
mask_BR = mask_Bd.reshape(-1,mask_Bd.shape[2])   # (T, ch)

phase_LR = phase_L.reshape(-1,phase_L.shape[2],phase_L.shape[3])   # (T, ch, freq)
mask_LR = mask_Ld.reshape(-1,mask_Ld.shape[2])   # (T, ch)

phase_MR = phase_M.reshape(-1,phase_M.shape[2],phase_H.shape[3])   # (T, ch, freq)
mask_MR = mask_Md.reshape(-1,mask_Md.shape[2])   # (T, ch)

phase_HR = phase_H.reshape(-1,phase_H.shape[2],phase_H.shape[3])   # (T, ch, freq)
mask_HR = mask_Hd.reshape(-1,mask_Hd.shape[2])   # (T, ch)

phase_BR.shape, mask_BR.shape, mask_BR.shape[0]/250/20


""" PICK CELL FOR THE SPIKING COUNT """

idx = 0 # index for cell and LFP 
cell = idx_cell_hpc[idx]
# lfp_ch = spk[2][idx_cell_hpc[idx]] # lfp_channel
# print('- cell:', cell, '| lfp ch:', lfp_ch)
print('- cell:', cell, ', total number of cells in CA1 is: ', idx_cell_hpc.size)

print('For the 4 epochs: spk length {}, phase length {}'.format(
spk[0][cell].shape, phase_BR.shape[0]*4 + phase_BR.shape[0]/2*4))

### 3.D Get spike activity for each epoch (20 min of activity), for a single cell

fs = 250 # final sampling frequency for speed, phase, lfp
L = phase_BR.shape[0] # length of one trimmed epoch, 20 min 

# baseline, starting at min 5
shift = 5*60*fs # 5 min shift at the beginning of the epoch, in numb of points
spk_epoch_B = spk[0][cell][shift:shift+L]
# Low dose, starting at min 35
shift = 35*60*fs # 5 min shift at the beginning of the epoch, in numb of points
spk_epoch_L = spk[0][cell][shift:shift+L]
# Mid dose, starting at min 65
shift = 65*60*fs # 5 min shift at the beginning of the epoch, in numb of points
spk_epoch_M = spk[0][cell][shift:shift+L]
# High dose, starting at min 95
shift = 95*60*fs # 5 min shift at the beginning of the epoch, in numb of points
spk_epoch_H = spk[0][cell][shift:shift+L]

spk_epoch_H.shape
print('SPIKE COUNT for cell, both low speed and high speed trials') 
print('Across epochs: ', np.sum(spk[0][cell]))
print('Baseline',np.sum(spk_epoch_B))
print('Low dose:',np.sum(spk_epoch_L))
print('Mid dose:',np.sum(spk_epoch_M))
print('High dose:',np.sum(spk_epoch_H))

#%%

start_time = time.time()

ch_start = 0

nch = phase_BR.shape[1]
for ch in range(ch_start,nch):
    
    hist_arr_B, bin_edges_B = freq_phase_map(phase_BR, mask_BR, spk_epoch_B, ch)
    hist_arr_L, bin_edges_L = freq_phase_map(phase_LR, mask_LR, spk_epoch_L, ch)
    hist_arr_M, bin_edges_M = freq_phase_map(phase_MR, mask_MR, spk_epoch_M, ch)
    hist_arr_H, bin_edges_H = freq_phase_map(phase_HR, mask_HR, spk_epoch_H, ch)
    
    # normalization factor for the PDF frequency-phase
    norm_B = np.sum(hist_arr_B, axis=1)[0]
    norm_L = np.sum(hist_arr_L, axis=1)[0]
    norm_M = np.sum(hist_arr_M, axis=1)[0]
    norm_H = np.sum(hist_arr_H, axis=1)[0]
    
    norm = np.max([np.sum(hist_arr_B, axis=1)[0], np.sum(hist_arr_L, axis=1)[0], np.sum(hist_arr_M, axis=1)[0], np.sum(hist_arr_H, axis=1)[0]])
    
    time_taken = time.time() - start_time
    print(f"Time of execution {time_taken} in sec")
    
    # plot_freq_phase_map_one_epoch(hist_arr_B, bin_edges_B, freq, 'baseline')
        
    # Create dictionaries to store normalization, bin values, histogram values
    norm_dict = {'B':norm_B, 'L':norm_L, 'M':norm_M, 'H':norm_H}   
    bin_dict = {'B': bin_edges_B, 'L': bin_edges_L, 'M': bin_edges_M, 'H': bin_edges_H}    
    hist_dict = {'B': hist_arr_B, 'L': hist_arr_L, 'M': hist_arr_M, 'H': hist_arr_H}
    
    # plot all the epochs together 
    plot_freq_phase_map_all_epochs(hist_dict, bin_dict, norm_dict, freq, rec, sess, cell, "oriens", "oriens", "SO", ch, save_flag = True)
    
    


#%%

variables_to_keep = ['spk']

# Now, we iterate over the list of global variables and delete those not in our list
for variable in list(globals()):
    if variable not in variables_to_keep:
        del globals()[variable]
    
    
