% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the PSD and Spectrogram for:
%   a) The whole HPC
%   b) Each stratum of the HPC
%   Input: pre-processed LFP (CSD) via "Ketamine_LFP_CSD_preprocessing.py"
%   Output: psd and spectrogram matlab files
%           plots/figure of psd and spectrograms for all HPC and each
%           stratum
%
% @ Gino Del Ferraro, Fenton Lab, NYU, Nov/Dec 2023

clear all; close all; clc;

%update path to your EEGLAB
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\');
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\sigprocfunc\'); 
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\guifunc\');
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\adminfunc\');
addpath('C:\Users\fentonlab\Desktop\luke\kfx\analysis\00-functions-tools\matlab-functions');
addpath('C:\Users\fentonlab\Desktop\luke\kfx\analysis\00-functions-tools\matlab-functions\save_and_load');
addpath('C:\Users\fentonlab\Desktop\luke\kfx\analysis\00-functions-tools\matlab-functions\plotting');
addpath('C:\Users\fentonlab\Desktop\luke\kfx\analysis\00-functions-tools\matlab-functions\other_tools');
addpath('C:\Users\fentonlab\Desktop\luke\kfx\analysis\00-functions-tools\matlab-tools');
addpath('C:\Users\fentonlab\Desktop\luke\kfx\analysis\00-functions-tools\matlab-tools\shadedErrorBars');

% INPUT PARAMETERS 
iSess = 8; % python session number 
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

% Session 0
% so = 1:5; % oriens layer  
% sp = 6:14; % pyramidal layer 
% rad = 15:17; % Radiatum 
% lm = 18:20; % Loc Mol 
% dup = 21:29; % Dentate upper part 
% dd = 30:38; % Dentate lower part

% Session 1
% so = 1:5;
% sp = 6:13;
% rad = 14:14;
% lm = 15:16;
% dup = 17:21;
% dd = 22:27;

% Session 2
% so = 1:5;
% sp = 6:11;
% rad = 12:17;
% lm = 18:22;
% dup = 23:25;
% dd = 26:28;

% Session 3
% so = 1:9;
% sp = 10:14;
% rad = 15:15;
% lm = 16:17;
% dup = 18:20;
% dd = 21:23;

% Session 8
so = 1:1; % there were no SO channels, so SO figures should be ignored.
sp = 2:7; % first channel was taken from SP and given to SO to avoid double-counting
rad = 8:18;
lm = 19:30;
dup = 31:36;
dd = 37:42;

% Session 16
% so = 1:5;
% sp = 6:11;
% rad = 12:12;
% lm = 13:13;
% dup = 14:23;
% dd = 24:33;

% WARNING! CODE BELOW HAS DUMMY CHANNEL VALUES, to run analysis across all HPC
% Session 17 (THESE CHANNEL RANGES ARE ARTIFICIAL, just to get all-HPC result figures)
% so = 1:6;    % from excel file
% sp = 7:13;   % from excel file
% rad = 14:16; % made up! ignore radiatum figures.
% lm = 17:19;  % made up! ignore loc mol figures.
% dup = 20:22; % made up! ignore dentate up figures.
% dd = 23:25;  % made up! ignore dendate down figures.

% Session 15 (THESE CHANNEL RANGES ARE ARTIFICIAL, just to get all-HPC result figures)
% so = 1:6;    % made up
% sp = 7:13;   % made up
% rad = 14:16; % made up
% lm = 17:19;  % made up
% dup = 20:22; % made up
% dd = 23:26;  % made up

% Session 10 (THESE CHANNEL RANGES ARE ARTIFICIAL, just to get all-HPC result figures)
% so = 1:1;    % from excel file
% sp = 2:12;   % from excel file
% rad = 13:15; % made up
% lm = 16:19;  % made up
% dup = 20:27; % made up
% dd = 28:35;  % made up

% Session 9 (THESE CHANNEL RANGES ARE ARTIFICIAL, just to get all-HPC result figures)
% so = 1:5;    % from excel file
% sp = 6:12;   % from excel file
% rad = 13:17; % made up
% lm = 18:24;  % made up
% dup = 25:36; % made up
% dd = 37:48;  % made up

% Session 7 (THESE CHANNEL RANGES ARE ARTIFICIAL, just to get all-HPC result figures)
% so = 1:21;   % from excel file (suspicious value)
% sp = 22:26;  % from excel file (suspicious value)
% rad = 27:27; % made up
% lm = 28:29;  % made up
% dup = 30:30; % made up
% dd = 31:31;  % made up

% Session 4 (THESE CHANNEL RANGES ARE ARTIFICIAL, just to get all-HPC result figures)
% so = 1:11;   % from excel file
% sp = 12:19;  % from excel file
% rad = 20:21; % made up
% lm = 22:22;  % made up
% dup = 23:23; % made up
% dd = 24:24;  % made up
% END CODE WITH DUMMY CHANNEL VALUES

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% directory and file names for the recording paths 
main_dir = 'C:\Users\fentonlab\Desktop\luke\LFPs\';
HPC_file = strcat(main_dir,'HPC_lfp_paths.mat');
PFC_file = strcat(main_dir,'PFC_lfp_paths.mat');

main_dir_HPC = 'C:\Users\fentonlab\Desktop\luke\LFPs\HPC';
main_dir_PFC = 'C:\Users\fentonlab\Desktop\luke\LFPs\PFC';
load(strcat(main_dir,'color_strata_v2.mat')); % load colors for plots

load(HPC_file) 
Paths_HPC = extract_paths(HPC_file_list); 
dir_rec = strcat(main_dir_HPC,Paths_HPC{sess}) % directory path  for the recording
dir_lfp = strcat(main_dir_HPC,Paths_HPC{sess},'\LFPs_and_masks') % directory path  for the recording

% load LFPs and masks 
[lfp, lfp_all, mask] = load_lfp_and_masks(dir_lfp,name_lfp);

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

% keyboard 

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


% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTROGRAMS      %%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters 
fs = 1250;
start = 1; % 30*fs;
ends = 60*fs;

% INPUT PARAMETERS for the Spectograms
min_start = 14; % starting minute for the computation of Spectrogram
min_end = 15;   % ending minute
% ----------------------------------


% Don't chance these values below 
spec_par.fs = 1250; % sampling rate 
spec_par.fk = [0 100]; % freq range
spec_par.tapers = [0.6 5]; % Time resolution, Freq. resoluzion
k = floor(2*spec_par.tapers(1)*spec_par.tapers(2) - 1); % number of tapers used
spec_par.dn = 0.05; % sliding step
% ----------------------------------

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
save_spectrograms(spec_rec, dir_rec, sprintf('spec_rec_HPC_lfp_CDS_reg_min_%d_%d.mat',min_start,min_end))

% load spectrograms whole HPC 
% spec_rec = load_spectrograms(dir_rec);


% %%%%%%%%%%%%%%%%%%%%%%%
% PLOT SPECTROGRAMS 
% %%%%%%%%%%%%%%%%%%%%%%%

% INPUT PARAMETERS for the Plots
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

plot_spectrograms_all_epochs_one_region(spec_rec, mask, step_t, step_f, minute, min_lab, 'RS Ket - Oriens', dir_rec, 'so', 1)
plot_spectrograms_all_epochs_one_region(spec_rec, mask, step_t, step_f, minute, min_lab, 'RS Ket - Pyramidale', dir_rec, 'sp', 1)
plot_spectrograms_all_epochs_one_region(spec_rec, mask, step_t, step_f, minute, min_lab, 'RS Ket - Radiatum', dir_rec, 'rad', 1)
plot_spectrograms_all_epochs_one_region(spec_rec, mask, step_t, step_f, minute, min_lab, 'RS Ket - LocMol', dir_rec, 'lm', 1)
plot_spectrograms_all_epochs_one_region(spec_rec, mask, step_t, step_f, minute, min_lab, 'RS Ket - Dentate Up', dir_rec, 'dup', 1)
plot_spectrograms_all_epochs_one_region(spec_rec, mask, step_t, step_f, minute, min_lab, 'RS Ket - Dentate Down', dir_rec, 'dd', 1)



