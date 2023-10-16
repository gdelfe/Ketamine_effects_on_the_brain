
%plot EEG from neuropixels
clear all; close all; clc;

%update path to your EEGLAB
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\');
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\sigprocfunc\'); 
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\guifunc\');
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\adminfunc\');
addpath('C:\Users\fentonlab\Desktop\Gino\chronux_2_11\');
addpath('C:\Users\fentonlab\Desktop\Gino\Gino_codes\');

main_dir = 'C:\Users\fentonlab\Desktop\Gino\LFPs\';
HPC_file = strcat(main_dir,'HPC_lfp_paths.mat');
PFC_file = strcat(main_dir,'PFC_lfp_paths.mat');

main_dir_HPC = 'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC';
main_dir_PFC = 'C:\Users\fentonlab\Desktop\Gino\LFPs\PFC';

load(HPC_file) 
Paths_HPC = extract_paths(HPC_file_list);

sess = 1;

BRAIN_reg_rec_dir = strcat(main_dir_HPC,Paths_HPC{sess})

display(['Loading LFP data ...'])
load(strcat(BRAIN_reg_rec_dir,'\lfp_B_epoch_low_speed.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\lfp_L_epoch_low_speed.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\lfp_M_epoch_low_speed.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\lfp_H_epoch_low_speed.mat')); % channel x minute x trial x lfp

% plot 1 min lfp to check data is loaded correctly 
figure;
plot(lfp_B{1,1}(1,:)); hold on 
plot(lfp_L{1,1}(1,:)); hold on 
plot(lfp_M{1,1}(1,:)); hold on 
plot(lfp_H{1,1}(1,:))

% aggregate trials 
[lfp_B_all, lfp_L_all, lfp_M_all, lfp_H_all] = aggregate_trials(lfp_B,lfp_L,lfp_M,lfp_H)

% compute spectrograms for each minute 
[spec_B, spec_L, spec_M, spec_H, f] = psd_across_min(lfp_B_all, lfp_L_all, lfp_M_all, lfp_H_all, 3, 100);

% plotting 
plot_psd_20_min(spec_B, spec_L, spec_M, spec_H, f, 'PSD all HPC - Stationary - RS Ketamine')
% plotting normalized 
plot_psd_20_min_normalize(spec_B, spec_L, spec_M, spec_H, f, 'PSD all HPC normalized - Stationary - RS Ketamine')


fs = 1250;
fk = [0 100]; % freq range
tapers = [0.4 5]; % Time resolution, Freq. resoluzion 
k = floor(2*tapers(1)*tapers(2) - 1) % number of tapers used 
dn = 0.005; % sliding step
fs = 1024; % sampling 
pad = 2;

[spec,f,ti] = compute_plot_spec(lfp_B_all{1}, tapers , fs, dn, fk, pad, 0.05,'min 1')













figure;
plot(f,log10(spec_B(min,:)),'LineWidth', 2); hold on
plot(f,log10(spec_L(min,:)),'LineWidth', 2); hold on
plot(f,log10(spec_M(min,:)),'LineWidth', 2); hold on
plot(f,log10(spec_H(min,:)),'LineWidth', 2); hold on

