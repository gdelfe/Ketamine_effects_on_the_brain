

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
ch = 1; min = 13; % channel and minute to look at 
W  = 3; % frequency resolution PSD
fk = 100; % max frequency PSD

% name files to load/save (depending on the re-referencing technique used)
name_lfp = '_CSD';
name_file_psd = 'psd_hpc_csd';
method = 'CSD';


% directory and file names for the recording paths 
main_dir = 'C:\Users\fentonlab\Desktop\Gino\LFPs\';
HPC_file = strcat(main_dir,'HPC_lfp_paths.mat');
PFC_file = strcat(main_dir,'PFC_lfp_paths.mat');

main_dir_HPC = 'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC';
main_dir_PFC = 'C:\Users\fentonlab\Desktop\Gino\LFPs\PFC';

load(HPC_file) 
Paths_HPC = extract_paths(HPC_file_list); 

dir_rec = strcat(main_dir_HPC,Paths_HPC{3}) % directory path  for the recording 
psd_1 = load_psd(dir_rec, name_file_psd);

dir_rec = strcat(main_dir_HPC,Paths_HPC{4}) % directory path  for the recording 
psd_2 = load_psd(dir_rec, name_file_psd);

plot_psd_20_min(psd_1.HPC, ['HPC ',method,' - PSD not normalized - RS Ketamine'],dir_rec,['HPC_',method],0,method)
plot_psd_20_min(psd_2.HPC, ['HPC ',method,' - PSD not normalized - RS Ketamine'],dir_rec,['HPC_',method],0,method)

plot_psd_20_min_1st_and_boost(psd_1.HPC, psd_2.HPC, '1st day vs boost', dir_rec, 'HPC', 1, method)
plot_psd_20_min_1st_and_boost_normalized(psd_1.HPC, psd_2.HPC, '1st day vs boost - normalized', dir_rec, 'HPC', 1, method)


