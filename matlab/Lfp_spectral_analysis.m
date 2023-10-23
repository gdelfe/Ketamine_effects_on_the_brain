
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
BRAIN_reg_rec_dir = strcat(main_dir_HPC,Paths_HPC{sess}) % directory path  for the recording 

% load LFP and masks 
load_lfp_and_masks(BRAIN_reg_rec_dir)


% plot 1 min lfp to check data is loaded correctly 
figure;
plot(lfp_B{1,1}(1,:)); hold on 
plot(lfp_L{1,1}(1,:)); hold on 
plot(lfp_M{1,1}(1,:)); hold on 
plot(lfp_H{1,1}(1,:))

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSD  %%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%

% aggregate trials 
[lfp_B_agg, lfp_L_agg, lfp_M_agg, lfp_H_agg] = aggregate_trials(lfp_B,lfp_L,lfp_M,lfp_H);

% compute PSD for each minute, for 20 min 
[psd_B, psd_L, psd_M, psd_H, f] = psd_across_min(lfp_B_agg, lfp_L_agg, lfp_M_agg, lfp_H_agg, 3, 100);
% save psd for 20 min 
save_psd(psd_B, psd_L, psd_M, psd_H, f,BRAIN_reg_rec_dir)

% plotting PSD
plot_psd_20_min(psd_B, psd_L, psd_M, psd_H, f, 'PSD all HPC - Stationary - RS Ketamine')
% plotting PSD normalized 
plot_psd_20_min_normalize(psd_B, psd_L, psd_M, psd_H, f, 'PSD all HPC normalized - Stationary - RS Ketamine')



% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTROGRAMS      %%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters 
start = 1 % 30*fs;
ends = 60*fs;
fs = 1250;
spec_par.fs = 1250;
spec_par.fk = [0 100]; % freq range
spec_par.tapers = [0.6 5]; % Time resolution, Freq. resoluzion
k = floor(2*spec_par.tapers(1)*spec_par.tapers(2) - 1) % number of tapers used
spec_par.dn = 0.05; % sliding step


% compute spectrograms whole HPC
[spec_rec, spec_tf] = compute_spectrograms_whole_rec(lfp_B_all,lfp_L_all,lfp_M_all,lfp_H_all,start,ends,spec_par);
% save spectrograms whole HPC
save_spectrograms(spec_rec, BRAIN_reg_rec_dir)

% load spectrograms whole HPC 
load_spectrograms(BRAIN_reg_rec_dir)

title_spec = sprintf('min %d',min);
step_t = 10;
step_f = 20;
epoch = 'high';

plot_spectrogram(X, spec, f, ti, fs, step_t, step_f, title_spec, epoch,[]);
plot_spectrogram(X, spec, f, ti, fs, step_t, step_f, title_spec, epoch, "zscore");


X = sq(lfp_B_all(min,start:ends,:))';
plot_20_min_spectrograms(spec_rec_H,X, spec_tf, fs, step_t, step_f, 1:10)


f_gamma = find(f>20 & f<50);












figure;
plot(f,log10(psd_B(min,:)),'LineWidth', 2); hold on
plot(f,log10(psd_L(min,:)),'LineWidth', 2); hold on
plot(f,log10(psd_M(min,:)),'LineWidth', 2); hold on
plot(f,log10(psd_H(min,:)),'LineWidth', 2); hold on

