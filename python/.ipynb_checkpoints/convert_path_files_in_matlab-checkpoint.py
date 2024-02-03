# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 12:47:40 2023

@author: Gino Del Ferraro, Oct 23

Convert pickle python files into matlab files.
The files contains the file path for each recording in the Ketamine experiment
"""

import pickle
from scipy.io import savemat
import os 

binFullPath = r'C:\Users\fentonlab\Desktop\Gino\LFPs'
HPC_path_file = os.path.join(r'C:\Users\fentonlab\Desktop\Gino\LFPs','HPC_lfp_paths.file')
PFC_path_file = os.path.join(r'C:\Users\fentonlab\Desktop\Gino\LFPs','PFC_lfp_paths.file')

# file names in matlab
HPC_path_file_MAT = os.path.join(r'C:\Users\fentonlab\Desktop\Gino\LFPs','HPC_lfp_paths.mat')
PFC_path_file_MAT = os.path.join(r'C:\Users\fentonlab\Desktop\Gino\LFPs','PFC_lfp_paths.mat')

 
with open(HPC_path_file,'rb') as f:
    HPC_file_list = pickle.load(f)

HPC_file_list = [str(path) for path in HPC_file_list]


# Convert the list to a dictionary before saving
data_dict = {'HPC_file_list': HPC_file_list}
savemat(HPC_path_file_MAT, data_dict)


with open(PFC_path_file,'rb') as f:
    PFC_file_list = pickle.load(f)

PFC_file_list = [str(path) for path in PFC_file_list]


# Convert the list to a dictionary before saving
data_dict = {'PFC_file_list': PFC_file_list}
savemat(PFC_path_file_MAT, data_dict)