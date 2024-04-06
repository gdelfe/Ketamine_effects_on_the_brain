%plot EEG from neuropixels
clear; close all; clc;

%update path to your EEGLAB
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\');
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\sigprocfunc\'); 
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\guifunc\');
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\adminfunc\');

%path to your dat file
%datFile = 'D:\Ketamine_Recording\2022-07-18_16-05-09\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.0\continuous.dat';
%datFile = 'D:\Ketamine_Recording\SPIKE_GLX_DATA\KSOUT_GLX\0718\0718_10min\0715_Gain1000_500_External_Repeat1_g0_imec1\0715_Gain1000_500_External_Repeat1_g0_t0.imec1.ap.bin';
datFile ='Z:\f\fentonlab\RAWDATA\NeuroPix\Ketamine\2022-10-06-03-55-00_M026_R-ketamine_mPFC_HPC_1_3_10_30mpk\2022-10-06-03-55-00_M026_R-ketamine_mPFC_HPC_1_3_10_30mpk_g0_imec1\2022-10-06-03-55-00_M026_R-ketamine_mPFC_HPC_1_3_10_30mpk_g0_t0.imec1.lf.bin';
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


eeg = single(eeg);

%select subset of channels (every 4th)
%eeg = eeg(1:chSkip:end,:);
eeg = eeg(1:end,:);

%reverse order (tip last channel)
eeg = eeg(end:-1:1,:);

%remove mean of each channel
for i = 1:size(eeg,1)
    eeg(i,:) = eeg(i,:) - mean(eeg(i,:));
end

%remove common mean
m = mean(eeg,1);
for i = 1:size(eeg,1)
    eeg(i,:) = eeg(i,:) - m;
end
   

%winlength = time to display
%spacing = amplitude
eegplot(eeg, 'srate', eegFS, 'winlength',1, 'spacing', 100); % winlength number :data display time (sec)0.5 ,  5, ,,,,
%saveas(fig,'07-27-2022-M015_SAL')

