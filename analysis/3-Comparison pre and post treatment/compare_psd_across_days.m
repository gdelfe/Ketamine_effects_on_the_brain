

%plot EEG from neuropixels
clear all; close all; clc;

%update path to your EEGLAB
% addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\');
% addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\sigprocfunc\'); 
% addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\guifunc\');
% addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\adminfunc\');
addpath('C:\Users\fentonlab\Desktop\luke\kfx\analysis\00-functions-tools\matlab-functions');
addpath('C:\Users\fentonlab\Desktop\luke\kfx\analysis\00-functions-tools\matlab-functions\save_and_load');
addpath('C:\Users\fentonlab\Desktop\luke\kfx\analysis\00-functions-tools\matlab-functions\plotting');
addpath('C:\Users\fentonlab\Desktop\luke\kfx\analysis\00-functions-tools\matlab-functions\other_tools');
addpath('C:\Users\fentonlab\Desktop\luke\kfx\analysis\00-functions-tools\matlab-tools');
addpath('C:\Users\fentonlab\Desktop\luke\kfx\analysis\00-functions-tools\matlab-tools\shadedErrorBars');

iSess = 2; % python session number PRE
sess_pre = iSess + 1; % session number 
iSess = 3; % python session number POST
sess_post = iSess + 1; % session number 

ch = 1; min = 13; % channel and minute to look at 
W  = 3; % frequency resolution PSD
fk = 100; % max frequency PSD

% name files to load/save (depending on the re-referencing technique used)
name_file_psd = 'psd_hpc_csd';
method = 'CSD';


% directory and file names for the recording paths 
main_dir = 'C:\Users\fentonlab\Desktop\luke\LFPs\';
HPC_file = strcat(main_dir,'HPC_lfp_paths.mat');
PFC_file = strcat(main_dir,'PFC_lfp_paths.mat');
load(strcat(main_dir,'color_strata.mat')); % load colors for plots

main_dir_HPC = 'C:\Users\fentonlab\Desktop\luke\LFPs\HPC';
main_dir_PFC = 'C:\Users\fentonlab\Desktop\luke\LFPs\PFC';

load(HPC_file) 
Paths_HPC = extract_paths(HPC_file_list); 

dir_rec = strcat(main_dir_HPC,Paths_HPC{sess_pre}) % directory path  for the recording 
psd_1 = load_psd(dir_rec, name_file_psd);

dir_rec = strcat(main_dir_HPC,Paths_HPC{sess_post}) % directory path  for the recording 
psd_2 = load_psd(dir_rec, name_file_psd);

plot_psd_20_min(psd_1.HPC, ['HPC ',method,' - PSD not normalized - RS Ketamine'],dir_rec,['HPC_',method],0,method, colors.HPC)
plot_psd_20_min(psd_2.HPC, ['HPC ',method,' - PSD not normalized - RS Ketamine'],dir_rec,['HPC_',method],0,method, colors.HPC)


plot_psd_20_min_pre_post(psd_1.HPC, psd_2.HPC, 'PRE and POST', dir_rec, 'HPC', 1)
plot_psd_20_min_pre_post_normalized(psd_1.HPC, psd_2.HPC, 'PRE and POST - normalized', dir_rec, 'HPC', 1)


