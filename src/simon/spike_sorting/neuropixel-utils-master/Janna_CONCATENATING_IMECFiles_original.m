%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%Wrapping AP and LF files together 
%
%Janna - April 2020
%
%NYU
%
% https://djoshea.github.io/neuropixel-utils/imec_dathttps://djoshea.github.io/neuropixel-utils/imec_dataset/aset/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Housekeepig
close all; clear all; clc;

addpath('neuropixel-utils-markSplitOnly\');

%Path
Output_folder = fullfile('E:\HD_ADN\20220805_JA179_lADN_red_h2_neurotar_freeV_fixedV_dark_g0\full_20220805_JA179_lADN_red_h2_neurotar\20220805_JA179_lADN_red_h2_neurotar_conc.imec.ap.bin');% give your file a name 

%% MAKING THE IMEC FILE
% channelMapFile = ('X:\Users\Janna\Matlab\Kilosort-main\configFiles\hStripe_bottom_row_kilosortChanMap.mat'); % this is for hstripe configuration , it matters here which one 
channelMapFile = ('D:\Angelaki_Lab\20220830_MC038_lADN_red_h2_neurotar_fixedV_dark_g0\20220830_MC038_lADN_red_h2_neurotar_fixedV_dark_g0_imec0\hStripe_bottom_second_kilosortChanMap.mat'); % this is for hstripe2 configuration , it matters here which one 

imec1 = Neuropixel.ImecDataset('E:\HD_ADN\Arena\20220805_JA179_lADN_red_h2_R_LD1_g0\20220805_JA179_lADN_red_h2_R_LD1_g0_imec0\20220805_JA179_lADN_red_h2_R_LD1_g0_t0.imec0.ap.bin', 'channelMap', channelMapFile);
imec2 = Neuropixel.ImecDataset('E:\HD_ADN\20220805_JA179_lADN_red_h2_neurotar_freeV_fixedV_dark_g0\20220805_JA179_lADN_red_h2_neurotar_freeV_fixedV_dark_g0_imec0\20220805_JA179_lADN_red_h2_neurotar_freeV_fixedV_dark_g0_t0.imec0.ap.bin', 'channelMap', channelMapFile);

%% CONCENATE DIFFERENT FILES TO SORT THEM TOGETHER 

imecList = {imec1, imec2};

tic
imecOut = Neuropixel.ImecDataset.writeConcatenatedFileMatchGains(imecList, Output_folder);
toc;