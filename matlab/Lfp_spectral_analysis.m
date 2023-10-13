
%plot EEG from neuropixels
clear all; close all; clc;

%update path to your EEGLAB
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\');
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\sigprocfunc\'); 
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\guifunc\');
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\adminfunc\');
addpath('C:\Users\fentonlab\Desktop\Gino\chronux_2_11\');
addpath('C:\Users\fentonlab\Desktop\Gino\Gino_codes\');

% main_dir = 'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC\2022-08-01_04-30-00_M015_RSK_mPFC_HPC_3_10_30mpk\';
main_dir = 'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC\2022-07-27-07-41-00_M015_SAL_PFC_HPC_0_0_0mpk\';

display(['Loading LFP data ...')
load(strcat(main_dir,'lfp_B_epoch_low_speed.mat')) % channel x minute x trial x lfp
load(strcat(main_dir,'lfp_L_epoch_low_speed.mat')) % channel x minute x trial x lfp
load(strcat(main_dir,'lfp_M_epoch_low_speed.mat')) % channel x minute x trial x lfp
load(strcat(main_dir,'lfp_H_epoch_low_speed.mat')) % channel x minute x trial x lfp

% plot 1 min lfp to check data is loaded correctly 
figure;
plot(lfp_B{1,1}(1,:)); hold on 
plot(lfp_L{1,1}(1,:)); hold on 
plot(lfp_M{1,1}(1,:)); hold on 
plot(lfp_H{1,1}(1,:))

% aggregate trials 
[lfp_B_all, lfp_L_all, lfp_M_all, lfp_H_all] = aggregate_trials(lfp_B,lfp_L,lfp_M,lfp_H)

% compute spectrograms for each minute 
[spec_B, spec_L, spec_M, spec_H, f] = spectrograms_across_min(lfp_B_all, lfp_L_all, lfp_M_all, lfp_H_all, 3, 100);

% plotting 
plot_psd_20_min(spec_B, spec_L, spec_M, spec_H, f, 'PSD all HPC - Stationary - RS Ketamine')
% plotting normalized 
plot_psd_20_min_normalize(spec_B, spec_L, spec_M, spec_H, f, 'PSD all HPC normalized - Stationary - RS Ketamine')




















figure;
plot(f,log10(spec_B(min,:)),'LineWidth', 2); hold on
plot(f,log10(spec_L(min,:)),'LineWidth', 2); hold on
plot(f,log10(spec_M(min,:)),'LineWidth', 2); hold on
plot(f,log10(spec_H(min,:)),'LineWidth', 2); hold on

