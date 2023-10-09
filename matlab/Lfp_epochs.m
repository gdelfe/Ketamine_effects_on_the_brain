
%plot EEG from neuropixels
clear all; close all; clc;

%update path to your EEGLAB
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\');
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\sigprocfunc\'); 
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\guifunc\');
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\adminfunc\');
addpath('C:\Users\fentonlab\Desktop\Gino\chronux_2_11\');

main_dir = 'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC\2022-08-01_04-30-00_M015_RSK_mPFC_HPC_3_10_30mpk\';

load(strcat(main_dir,'lfp_B_epoch_low_speed.mat'))