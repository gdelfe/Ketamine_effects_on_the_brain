# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 16:11:37 2024

@author: Gino Del Ferraro

"""
import sys
import os 

dir_current = os.path.dirname(__file__)
dir_utils = os.path.join(dir_current, '..', 'utils-tools')
sys.path.append(dir_utils)

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


binFullPath = r'C:\Users\fentonlab\Desktop\Gino\LFPs'
main_dir = r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC'
HPC_path_file = os.path.join(r'C:\Users\fentonlab\Desktop\Gino\LFPs','HPC_lfp_paths.file')
PFC_path_file = os.path.join(r'C:\Users\fentonlab\Desktop\Gino\LFPs','PFC_lfp_paths.file')

sess = 2

# ====== Load rec path 
rec = load_rec_path(binFullPath,HPC_path_file,PFC_path_file,"HPC",sess)
idx_cell_HPC = load_list_channel_CA1(rec, sess, "HPC")

idx = 0
cell = idx_cell_hpc[idx]
#%%

# frequency-phase map, list of 4 elements, one for each layer, each layer has 4 epochs
# call it via pdf_sm[0]["B"]

pdf_sm = load_matrices_freq_phase_4_by_4(rec, sess, "HPC", cell, idx) 
