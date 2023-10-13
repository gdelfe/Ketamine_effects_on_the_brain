
%plot EEG from neuropixels
clear all; close all; clc;

%update path to your EEGLAB
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\');
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\sigprocfunc\'); 
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\guifunc\');
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\adminfunc\');
addpath('C:\Users\fentonlab\Desktop\Gino\chronux_2_11\');

main_dir = 'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC\2022-08-01_04-30-00_M015_RSK_mPFC_HPC_3_10_30mpk\';

load(strcat(main_dir,'lfp_B_epoch_low_speed.mat')) % channel x minute x trial x lfp
load(strcat(main_dir,'lfp_L_epoch_low_speed.mat')) % channel x minute x trial x lfp
load(strcat(main_dir,'lfp_M_epoch_low_speed.mat')) % channel x minute x trial x lfp
load(strcat(main_dir,'lfp_H_epoch_low_speed.mat')) % channel x minute x trial x lfp


lfp_B_all = cell(size(lfp_B,2),1); % cell min x
lfp_L_all = cell(size(lfp_L,2),1); % cell min x
lfp_M_all = cell(size(lfp_M,2),1); % cell min x
lfp_H_all = cell(size(lfp_H,2),1); % cell min x

% aggregate trials for different channels together, minute by minute 
for min = 1:size(lfp_B,2) % for all the min 
    for ch = 1:size(lfp_B,1) % for all the channels 
        lfp_B_all{min} = cat(1,lfp_B_all{min},lfp_B{ch,min}); % stack along 1st dimension all the channels for the same minute, final shape: min x trial x lfp values (1 min)
        lfp_L_all{min} = cat(1,lfp_L_all{min},lfp_L{ch,min}); % stack along 1st dimension all the channels for the same minute
        lfp_M_all{min} = cat(1,lfp_M_all{min},lfp_M{ch,min}); % stack along 1st dimension all the channels for the same minute
        lfp_H_all{min} = cat(1,lfp_H_all{min},lfp_H{ch,min}); % stack along 1st dimension all the channels for the same minute
    end
end 












