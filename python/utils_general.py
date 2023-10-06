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
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import pandas as pd
import os
import h5py
import warnings
import copy
from scipy.stats import zscore
from urllib.request import urlopen
plt.rcParams['pdf.fonttype'] = 42

#import mne
import sys
import pickle


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
            
        for ch, path in enumerate(HPC_file_list):
            rec[ch]['HPC'] =  str(path)
    
    elif brain_reg == 'PFC':
        # load PFC file names and store them in rec 
        with open(PFC_path_file,'rb') as f:
            PFC_file_list = pickle.load(f)
            
        for ch, path in enumerate(PFC_file_list):
            rec[ch]['PFC'] =  str(path)
        
    else:
        sys.exit('Brain regio is not HPC or PFC -- exit')
        
        
    # select path for specific session recording 
    binFullPath = rec[sess][brain_reg] 
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
        print('Gain factor multiplication should be done channel by channel, this is just a bit slower')
        
        
    #### Read LFP start/end time, HPC channels - Choose whether loading HPC or PFC here
        
    pd.set_option('display.max_colwidth',100)
    in_file = os.path.join(r'C:\Users\fentonlab\Desktop\Gino\LFPs\\', "lfp_aln_idx.csv")
    Lfp_aln = pd.read_csv(in_file)
    

    print(Lfp_aln.head())
    
    
    # =============================================================================
    #  Trim LFP 
    # =============================================================================
    #### Pull Lfp 'start' and 'end' for the specific recording 
    
    start = Lfp_aln["Lfp start"][sess] # starting of the trial in Lfp data points
    end = Lfp_aln["Lfp end"][sess] # end of trial
    # ch_start = Lfp_aln["HPC ch start"][sess]  # Hpc first channel
    # ch_end = Lfp_aln["HPC ch end"][sess]  # Hpc last channel
    
    # 0-120 CA1
    # 0-86 ripples 
    # 126-244 dentate
    
    
    ch_start = 0
    ch_end = 120
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

    return Lfp, speed 



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
    Lfp_B.shape, speed_B.shape
    
    return Lfp_B, Lfp_L, Lfp_M, Lfp_H, speed_B, speed_L, speed_M, speed_H

# =============================================================================


def select_1min_data(current_min, N=2500):
    
    #### 3. Select 1 min of data for each Epoch, for both Lfp and speed
    ##### Select 1 every M channels for the Lfp

    # parameters
    current_min = 0 # minute to look at

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


    Lfp_B_min.shape, speed_B_min.shape, speed_B_min.size/N
    
    



