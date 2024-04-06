# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 16:11:37 2024

@author: Gino Del Ferraro

"""
import sys
import os 

dir_current = os.path.dirname(__file__)
dir_utils = os.path.abspath(os.path.join(dir_current, '..', 'utils-tools'))

print(dir_utils)
sys.path.append(dir_utils)


import numpy as np
import pickle
from utilities_ketamine_analysis_v8 import *
from utils_signal_processing import *
from utils_plotting import *
from utils_general import *
from utils_phase_frequency import *
import matplotlib.gridspec as gridspec
import time 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

binFullPath = r'C:\Users\fentonlab\Desktop\Gino\LFPs'
main_dir = r'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC'
HPC_path_file = os.path.join(r'C:\Users\fentonlab\Desktop\Gino\LFPs','HPC_lfp_paths.file')
PFC_path_file = os.path.join(r'C:\Users\fentonlab\Desktop\Gino\LFPs','PFC_lfp_paths.file')

sess = 2

# ====== Load rec path 
rec = load_rec_path(binFullPath,HPC_path_file,PFC_path_file,"HPC",sess)
idx_cell_HPC = load_list_channel_CA1(rec, sess, "HPC")
freq = np.logspace(0.32, 2.48, num=70, endpoint=True)

# frequency range/bands
theta = np.where((freq>=2) & (freq <=15)) # [2, 15] Hz
low_gamma = np.where((freq>=30) & (freq <=60)) # [30, 60] Hz
high_gamma = np.where((freq>60) & (freq <=100)) # # (60, 100] Hz

"""
Compute Mean Abs Error for freq-phase map wrt to Baseline-Low dose, Baseline-Mid dose, Baseline-High dose
in the theta, low gamma, and high gamma frequency range 
"""
mae_sess = []
for idx in range(0,len(idx_cell_HPC)):
    cell = idx_cell_HPC[idx]
    pdf_sm = load_matrices_freq_phase_4_by_4(rec, sess, "HPC", cell, idx) 
    mae_sess.append( mean_abs_error_freq_phase_maps(pdf_sm,theta,low_gamma,high_gamma) )
    

# vector = create_array_for_3D_plots(mae_sess,idx_cell_HPC,0,'theta')
# plot_3D_mean_abs_error_all_cells(vector)

vector_freq = create_array_for_3D_plots_axis_freq_band(mae_sess,idx_cell_HPC,0,'BL')
plot_3D_mean_abs_error_all_cells(vector_freq,'Baseline - Low dose')

vector_freq = create_array_for_3D_plots_axis_freq_band(mae_sess,idx_cell_HPC,0,'BM')
plot_3D_mean_abs_error_all_cells(vector_freq,'Baseline - Mid dose')

vector_freq = create_array_for_3D_plots_axis_freq_band(mae_sess,idx_cell_HPC,0,'BH')
plot_3D_mean_abs_error_all_cells(vector_freq,'Baseline - High dose')




# plot_3D_mean_abs_error_all_cells_animated(vector, 'animatin')



