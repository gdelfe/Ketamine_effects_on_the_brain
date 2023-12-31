
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
ch = 1; minute = 13; % channel and minute to look at 
W  = 3; % frequency resolution PSD
fk = 100; % max frequency PSD

% name files to load/save (depending on the re-referencing technique used)
name_lfp = '_CSD';
name_file_psd = 'psd_hpc_csd';
method = 'CSD';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channels for HPC subareas -- Copy/Pase the values on the excel file
% Session 2
so = 1:5;  % oriens layer  
sp = 6:11; % pyramidal layer 
rad = 12:17; % Radiatum 
lm = 18:22; % Loc Mol 
dup = 23:25; % Dentate upper part 
dd = 26:28; % Dentate lower part
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channels for HPC subareas -- Copy/Pase the values on the excel file
% Session 3
% CA1 = 1:9;  % oriens 
% ripple = 10:14; % pyramidal layer 
% rad = 15:15; % Radiatum 
% lm = 16:17; % Loc Mol 
% dup = 18:20 % Dentate upper part 
% dd = 21:23; % Dentate lower part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% directory and file names for the recording paths 
main_dir = 'C:\Users\fentonlab\Desktop\Gino\LFPs\';
HPC_file = strcat(main_dir,'HPC_lfp_paths.mat');
PFC_file = strcat(main_dir,'PFC_lfp_paths.mat');

main_dir_HPC = 'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC';
main_dir_PFC = 'C:\Users\fentonlab\Desktop\Gino\LFPs\PFC';
load(strcat(main_dir,'color_strata.mat')); % load colors for plots

load(HPC_file) 
Paths_HPC = extract_paths(HPC_file_list); 
dir_rec = strcat(main_dir_HPC,Paths_HPC{sess}) % directory path  for the recording 

% load LFPs and masks 
[lfp, lfp_all, mask] = load_lfp_and_masks(dir_rec,name_lfp);

nch = size(lfp_all.B,3);
% plot 1 min lfp to check data is loaded correctly 
figure;
plot(lfp.B{ch,minute}(1,:)); hold on 
plot(lfp.L{ch,minute}(1,:)); hold on 
plot(lfp.M{ch,minute}(1,:)); hold on 
plot(lfp.H{ch,minute}(1,:))

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSD  %%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%

% aggregate trials 
[lfp_agg.HPC] = aggregate_trials_ch_range(lfp,1:nch); % All HPC
[lfp_agg.so] = aggregate_trials_ch_range(lfp,so);
[lfp_agg.sp] = aggregate_trials_ch_range(lfp,sp);
[lfp_agg.rad] = aggregate_trials_ch_range(lfp,rad);
[lfp_agg.lm] = aggregate_trials_ch_range(lfp,lm);
[lfp_agg.dup] = aggregate_trials_ch_range(lfp,dup);
[lfp_agg.dd] = aggregate_trials_ch_range(lfp,dd);

% compute PSD for each minute, for 20 min 
[psd.HPC] = psd_across_min(lfp_agg.HPC, W, fk); % All HPC
[psd.so] = psd_across_min(lfp_agg.so, W, fk);
[psd.sp] = psd_across_min(lfp_agg.sp, W, fk);
[psd.rad] = psd_across_min(lfp_agg.rad, W, fk);
[psd.lm] = psd_across_min(lfp_agg.lm, W, fk);
[psd.dup] = psd_across_min(lfp_agg.dup, W, fk);
[psd.dd] = psd_across_min(lfp_agg.dd, W, fk);

% save psd for 20 min 
save_psd(psd, dir_rec, name_file_psd);
% psd = load_psd(dir_rec, name_file_psd);

keyboard 

% plotting PSD
% plot_psd_20_min(psd, 'PSD all HPC - Stationary - RS Ketamine', dir_rec,1)
% plotting PSD normalized 
% All HPC 
plot_psd_20_min(psd.HPC, ['HPC ',method,' - PSD not normalized - RS Ketamine'],dir_rec,['HPC_',method],1,method, colors.HPC)
plot_psd_20_min_normalize(psd.HPC, ['HPC ',method,' - PSD normalized - RS Ketamine'],dir_rec,['HPC_',method],1,method, colors.HPC)
% By region 
plot_psd_20_min_normalize(psd.so, ['S. Oriens ',method,' - PSD normalized - RS Ketamine'],dir_rec,['SO_',method],1,method, colors.so)
plot_psd_20_min_normalize(psd.sp, ['S. Pyramidale ',method ,' - PSD normalized  RS Ketamine'], dir_rec, ['SP_',method],1,method,  colors.sp)
plot_psd_20_min_normalize(psd.rad, ['Radiatum ',method,' - PSD normalized - RS Ketamine'], dir_rec, ['Radiatum_',method],1,method, colors.rad)
plot_psd_20_min_normalize(psd.lm, ['LocMol ',method,' - PSD normalized - RS Ketamine'], dir_rec, ['LocMol_',method],1,method, colors.lm)
plot_psd_20_min_normalize(psd.dup, ['Dentate Up ',method,' - PSD normalized  - RS Ketamine'], dir_rec,['DentUp_',method], 1,method, colors.dup)
plot_psd_20_min_normalize(psd.dd, ['Dentate Down ',method,' - PSD normalized - RS Ketamine'], dir_rec, ['DentDown_',method],1,method, colors.dd)

% plot 1 min PSD for each region, each epoch 
minute = 10;
plot_psd_all_HPC_region_all_epochs_compact(psd, dir_rec, minute, 1, 'CSD')
plot_psd_all_HPC_region_all_epochs_normalized_compact(psd, dir_rec, minute, 1, 'CSD')

keyboard 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTROGRAMS      %%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters 
fs = 1250;
start = 1; % 30*fs;
ends = 60*fs;
min_start = 14;
min_end = 15;

spec_par.fs = 1250;
spec_par.fk = [0 100]; % freq range
spec_par.tapers = [0.6 5]; % Time resolution, Freq. resoluzion
k = floor(2*spec_par.tapers(1)*spec_par.tapers(2) - 1) % number of tapers used
spec_par.dn = 0.05; % sliding step;

% compute spectrograms whole HPC
spec_rec = {};
[spec_rec] = compute_spectrograms_whole_rec(spec_rec, lfp_all, 'HPC', fs, min_start, min_end, start, ends, spec_par, 1:nch); % all HPC 
[spec_rec] = compute_spectrograms_whole_rec(spec_rec, lfp_all, 'so', fs, min_start, min_end,start, ends, spec_par, so);
[spec_rec] = compute_spectrograms_whole_rec(spec_rec, lfp_all, 'sp',fs, min_start, min_end, start, ends, spec_par, sp);
[spec_rec] = compute_spectrograms_whole_rec(spec_rec, lfp_all,'rad', fs, min_start, min_end, start, ends, spec_par, rad);
[spec_rec] = compute_spectrograms_whole_rec(spec_rec, lfp_all,'lm', fs, min_start, min_end, start, ends, spec_par, lm);
[spec_rec] = compute_spectrograms_whole_rec(spec_rec, lfp_all, 'dup', fs, min_start, min_end, start, ends, spec_par, dup);
[spec_rec] = compute_spectrograms_whole_rec(spec_rec, lfp_all, 'dd', fs, min_start, min_end, start, ends, spec_par, dd);

% save spectrograms whole HPC
save_spectrograms(spec_rec, dir_rec, '\spec_rec_HPC_lfp_CDS_reg_min_14_15.mat')

% load spectrograms whole HPC 
% spec_rec = load_spectrograms(dir_rec);


% %%%%%%%%%%%%%%%%%%%%%%%
% PLOT SPECTROGRAMS 
% %%%%%%%%%%%%%%%%%%%%%%%

step_t = 10;
step_f = 20;
minute = 2;       % This is the min in spec, depends on what min you are looking at 
min_lab = 15;  % this depends on min_lab above

% plot single spectrograms:
title_spec = sprintf('min %d',minute);
epoch = 'high';
% no-normalization
plot_spectrogram(sq(spec_rec.H), 'HPC', spec_rec, minute, step_t, step_f, title_spec, epoch,[], dir_rec, 1);
% zscore normalization 
plot_spectrogram(sq(spec_rec.H),'HPC', spec_rec, minute, step_t, step_f, title_spec, epoch, "zscore",dir_rec, 1);

% plot multiple spectrograms, same epoch 
range = 1:10;
% plot_20_min_spectrograms(spec_rec.H.HPC, spec_rec, step_t, step_f, range,'high','RS Ket',dir_rec,1)

% plot 1 min spectrogram, all epochs 
plot_spectrograms_all_epochs(spec_rec,'HPC', mask, step_t, step_f, minute, min_lab, 'RS Ket',dir_rec,1)
% plot_spectrograms_all_epochs_and_gamma(spec_rec, mask, step_t, step_f, min,'RS Ket',dir_rec,1)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SPECTROGRAMS FOR SUBREGIONS HPC
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLOT SPECTROGRAM FOR EACH EPOCH, EACH SUBREGION 
plot_spectrograms_all_regions(spec_rec, 'B', mask, step_t, step_f, minute, min_lab, 'RS Ket - BASELINE', dir_rec, 1)
plot_spectrograms_all_regions(spec_rec, 'L', mask, step_t, step_f, minute, min_lab, 'RS Ket - LOW DOSE', dir_rec, 1)
plot_spectrograms_all_regions(spec_rec, 'M', mask, step_t, step_f, minute, min_lab, 'RS Ket - MID DOSE', dir_rec, 1)
plot_spectrograms_all_regions(spec_rec, 'H', mask, step_t, step_f, minute, min_lab, 'RS Ket - HIGH DOSE', dir_rec, 1)

plot_spectrograms_all_epochs_one_region(spec_rec, mask, step_t, step_f, minute, min_lab, 'RS Ket - HPC', dir_rec, 'HPC', 1)

plot_spectrograms_all_epochs_one_region(spec_rec, mask, step_t, step_f, minute, min_lab, 'RS Ket - CA1', dir_rec, 'CA1', 1)
plot_spectrograms_all_epochs_one_region(spec_rec, mask, step_t, step_f, minute, min_lab, 'RS Ket - Ripple', dir_rec, 'ripple', 1)
plot_spectrograms_all_epochs_one_region(spec_rec, mask, step_t, step_f, minute, min_lab, 'RS Ket - Radiatum', dir_rec, 'rad', 1)
plot_spectrograms_all_epochs_one_region(spec_rec, mask, step_t, step_f, minute, min_lab, 'RS Ket - LocMol', dir_rec, 'lm', 1)
plot_spectrograms_all_epochs_one_region(spec_rec, mask, step_t, step_f, minute, min_lab, 'RS Ket - Dentate Up', dir_rec, 'dup', 1)
plot_spectrograms_all_epochs_one_region(spec_rec, mask, step_t, step_f, minute, min_lab, 'RS Ket - Dentate Down', dir_rec, 'dd', 1)

plot_zscored_psd_from_spectrogram(spec_rec,dir_rec,minute,min_lab,1)
% plot_psd_from_spectrogram(spec_rec,dir_rec,min,min_lab,1)










figure;
plot(f,log10(psd_B(minute,:)),'LineWidth', 2); hold on
plot(f,log10(psd_L(minute,:)),'LineWidth', 2); hold on
plot(f,log10(psd_M(minute,:)),'LineWidth', 2); hold on
plot(f,log10(psd_H(minute,:)),'LineWidth', 2); hold on

