
%plot EEG from neuropixels
clear all; close all; clc;

%update path to your EEGLAB
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\');
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\sigprocfunc\'); 
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\guifunc\');
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\adminfunc\');
addpath('C:\Users\fentonlab\Desktop\Gino\chronux_2_11\');

main_dir = 'Z:\f\fentonlab\RAWDATA\NeuroPix\Ketamine\';
file_path = '2022-08-01_04-30-00_M015_RSK_mPFC_HPC_3_10_30mpk\2022-08-01_04-30-00_M015_RSK_mPFC_HPC_3_10_30mpk_g0_imec0\2022-08-01_04-30-00_M015_RSK_mPFC_HPC_3_10_30mpk_g0_t0.imec0.lf.bin';
datFile = strcat(main_dir,file_path);
%path to your dat file
%datFile = 'D:\Ketamine_Recording\2022-07-18_16-05-09\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.0\continuous.dat';
%datFile = 'D:\Ketamine_Recording\SPIKE_GLX_DATA\KSOUT_GLX\0718\0718_10min\0715_Gain1000_500_External_Repeat1_g0_imec1\0715_Gain1000_500_External_Repeat1_g0_t0.imec1.ap.bin';
% datFile ='Z:\f\fentonlab\RAWDATA\NeuroPix\Ketamine\2022-10-06-03-55-00_M026_R-ketamine_mPFC_HPC_1_3_10_30mpk\2022-10-06-03-55-00_M026_R-ketamine_mPFC_HPC_1_3_10_30mpk_g0_imec1\2022-10-06-03-55-00_M026_R-ketamine_mPFC_HPC_1_3_10_30mpk_g0_t0.imec1.lf.bin';
files = dir(datFile);

%eegFS = 30000;
eegFS = 2500; %downsampled FS with LFP %%AP 30000 Sample 30000*10sec

Nchannels = 385;

%chSkip = 4;

%Nbytes = files(1).bytes;
%Nsamples = Nbytes/Nchannels/2; %% for whole file
%Nsamples =eegFS*1; %%(10sec)   %% for short time checking
Nsamples =eegFS*10;

%open file
fid = fopen(datFile, 'r');
%read data
eeg = fread(fid, [Nchannels Nsamples], '*int16');
fclose(fid); %close file

lfp = double(eeg(:,1:Nsamples));

figure;
plot(lfp(35,:))


% %%%%%%%%%%%%%%%%
% PSD 
% %%%%%%%%%%%%%%%%

% parameters
fs = 2500; % sampling frequency
W = 5;  % freq resolution 
N = size(lfp,2)/fs; % Time length in sec
K = 2*N*W -1 % number of tapers 
fk = 200; % max freq

tapers = [N W];
pad = 2;


[spec, f] = dmtspec(lfp, tapers, fs, fk, pad, 0.05, 1)

figure;
semilogy(f,spec,'LineWidth', 2); hold on
xlim([0 20])







