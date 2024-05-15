
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
dir_rec = strcat(main_dir_HPC,Paths_HPC{sess}) % directory path  for the recording 

% load LFP and masks 
[lfp, lfp_all, mask] = load_lfp_and_masks(dir_rec);


% plot 1 min lfp to check data is loaded correctly 
ch = 1; min = 13;
figure;
plot(lfp.B{ch,min}(1,:)); hold on 
plot(lfp.L{ch,min}(1,:)); hold on 
plot(lfp.M{ch,min}(1,:)); hold on 
plot(lfp.H{ch,min}(1,:))

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSD  %%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%

% aggregate trials 
[lfp_agg] = aggregate_trials(lfp);

% compute PSD for each minute, for 20 min 
[psd] = psd_across_min(lfp_agg, 3, 100);

% save psd for 20 min 
save_psd(psd, dir_rec,'psd');
psd = load_psd(dir_rec, 'psd');


% plotting PSD
plot_psd_20_min(psd, 'All HPC - PSD - RS Ketamine', dir_rec,1)
% plotting PSD normalized 
plot_psd_20_min_normalize(psd, 'All HPC - PSD normalized - RS Ketamine',dir_rec,'all_HPC',1)


% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTROGRAMS      %%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters 
fs = 1250;
start = 1; % 30*fs;
ends = 60*fs;

spec_par.fs = 1250;
spec_par.fk = [0 100]; % freq range
spec_par.tapers = [0.6 5]; % Time resolution, Freq. resoluzion
k = floor(2*spec_par.tapers(1)*spec_par.tapers(2) - 1) % number of tapers used
spec_par.dn = 0.05; % sliding step;

% compute spectrograms whole HPC
[spec_rec] = compute_spectrograms_whole_rec(lfp_all, fs, start, ends, spec_par);

% save spectrograms whole HPC
save_spectrograms(spec_rec, dir_rec)

% load spectrograms whole HPC 
spec_rec = load_spectrograms(dir_rec);


% %%%%%%%%%%%%%%%%%%%%%%%
% PLOT SPECTROGRAMS 
% %%%%%%%%%%%%%%%%%%%%%%%

step_t = 10;
step_f = 20;

% plot single spectrograms:
title_spec = sprintf('min %d',min);
epoch = 'high';
% no-normalization
plot_spectrogram(sq(spec_rec.H(:,:,min)), spec_rec, step_t, step_f, title_spec, epoch,[],dir_rec, 1);
% zscore normalization 
plot_spectrogram(sq(spec_rec.H(:,:,min)), spec_rec, step_t, step_f, title_spec, epoch, "zscore",dir_rec, 1);

% plot multiple spectrograms, same epoch 
range = 10:16;
plot_20_min_spectrograms(spec_rec.M, spec_rec, step_t, step_f, range,'MID','RS Ket',dir_rec,1)

% plot 1 min spectrogram, all epochs 
hpc_area = 'all_HPC';
min = 15;
plot_spectrograms_all_epochs(spec_rec, mask, step_t, step_f, min,'RS Ket',dir_rec, hpc_area,1)
plot_spectrograms_all_epochs_and_gamma(spec_rec, mask, step_t, step_f, min,'RS Ket',dir_rec, hpc_area,1)





