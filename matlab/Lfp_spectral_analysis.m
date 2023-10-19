
%plot EEG from neuropixels
clear all; close all; clc;

%update path to your EEGLAB
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\');
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\sigprocfunc\'); 
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\guifunc\');
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\adminfunc\');
addpath('C:\Users\fentonlab\Desktop\Gino\chronux_2_11\');
addpath('C:\Users\fentonlab\Desktop\Gino\Gino_codes\');

iSess = 2; % python session number 
sess = iSess + 1; % session number 

% directory and file names for the recording paths 
main_dir = 'C:\Users\fentonlab\Desktop\Gino\LFPs\';
HPC_file = strcat(main_dir,'HPC_lfp_paths.mat');
PFC_file = strcat(main_dir,'PFC_lfp_paths.mat');

main_dir_HPC = 'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC';
main_dir_PFC = 'C:\Users\fentonlab\Desktop\Gino\LFPs\PFC';

load(HPC_file) 
Paths_HPC = extract_paths(HPC_file_list);
BRAIN_reg_rec_dir = strcat(main_dir_HPC,Paths_HPC{sess})

display(['Loading LFP data ...'])
%
load(strcat(BRAIN_reg_rec_dir,'\lfp_B_epoch_low_speed.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\lfp_L_epoch_low_speed.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\lfp_M_epoch_low_speed.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\lfp_H_epoch_low_speed.mat')); % channel x minute x trial x lfp

% lfp_X_all type of variables 
load(strcat(BRAIN_reg_rec_dir,'\lfp_B_epoch_all_trials.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\lfp_L_epoch_all_trials.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\lfp_M_epoch_all_trials.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\lfp_H_epoch_all_trials.mat')); % channel x minute x trial x lfp

% low speed masks
load(strcat(BRAIN_reg_rec_dir,'\mask_low_B.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\mask_low_L.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\mask_low_M.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\mask_low_H.mat')); % channel x minute x trial x lfp


% high speed masks
load(strcat(BRAIN_reg_rec_dir,'\mask_high_B.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\mask_high_L.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\mask_high_M.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\mask_high_H.mat')); % channel x minute x trial x lfp


% plot 1 min lfp to check data is loaded correctly 
figure;
plot(lfp_B{1,1}(1,:)); hold on 
plot(lfp_L{1,1}(1,:)); hold on 
plot(lfp_M{1,1}(1,:)); hold on 
plot(lfp_H{1,1}(1,:))

% aggregate trials 
[lfp_B_agg, lfp_L_agg, lfp_M_agg, lfp_H_agg] = aggregate_trials(lfp_B,lfp_L,lfp_M,lfp_H)

% compute PSD for each minute 
[spec_B, spec_L, spec_M, spec_H, f] = psd_across_min(lfp_B_agg, lfp_L_agg, lfp_M_agg, lfp_H_agg, 3, 100);

% plotting PSD
plot_psd_20_min(spec_B, spec_L, spec_M, spec_H, f, 'PSD all HPC - Stationary - RS Ketamine')
% plotting PSD normalized 
plot_psd_20_min_normalize(spec_B, spec_L, spec_M, spec_H, f, 'PSD all HPC normalized - Stationary - RS Ketamine')




fs = 1250;
start = 1 % 30*fs;
ends = 60*fs;
fk = [0 100]; % freq range
tapers = [0.6 5]; % Time resolution, Freq. resoluzion
k = floor(2*tapers(1)*tapers(2) - 1) % number of tapers used
dn = 0.05; % sliding step
pad = 2;

min = 15;
X = sq(lfp_B_all(min,start:ends,:))';
title_spec = sprintf('min %d',min);
step_t = 10;
step_f = 10;
epoch = 'high';

[spec, f, ti] = tfspec(X, tapers, fs, dn, fk, pad, 0.05,1);
plot_spectrogram(X, spec, f, ti, fs, step_t, step_f, title_spec, epoch,[]);
plot_spectrogram(X, spec, f, ti, fs, step_t, step_f, title_spec, epoch, "zscore");


















figure;
plot(f,log10(spec_B(min,:)),'LineWidth', 2); hold on
plot(f,log10(spec_L(min,:)),'LineWidth', 2); hold on
plot(f,log10(spec_M(min,:)),'LineWidth', 2); hold on
plot(f,log10(spec_H(min,:)),'LineWidth', 2); hold on

