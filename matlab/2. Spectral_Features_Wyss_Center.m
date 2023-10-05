
% Code number 2
%
% Machine Learning Scientist task
% 
% This code generates the spectrograms to be used in the classification
% problem. It also plots the avg and single trials spectrograms for a few 
% trials (without saving them)
%
% @ Gino Del Ferraro, NYU, June 2023


addpath('tools/')

clear all; close all

dir_data = '../Data/' 

% EVENTS BP + notch MNE

sess = 2; % Change this value to work on each one of the different sessions

% Load 2.5 sec time series for R, M1, M2
load(strcat(dir_data,sprintf('R_F_split_ord_sess_%d.mat',sess))) % Rest events 
load(strcat(dir_data,sprintf('M1_F_split_ord_sess_%d.mat',sess))) % M1 events
load(strcat(dir_data,sprintf('M2_F_split_ord_sess_%d.mat',sess))) % M2 events 

% (channel, trial, time)
RFch = RF;
M1ch = M1F;
M2ch = M2F; 

size(RFch)
size(M1ch)
size(M2ch)


% stack trials of different channels together for PSD calculation 
RF = reshape(RF,[],size(RF,3));
M1F = reshape(M1F,[],size(M1F,3));
M2F = reshape(M2F,[],size(M2F,3));


% %%%%%%%%%%%%%%%%%%%%%%%
% Power Spectral Density
% %%%%%%%%%%%%%%%%%%%%%%

W = 3; % Frequency resolution for PSD 
N = 2560; % Lenght of time series
sampling = 1024; % sampling 
fk = 100; % max frequency

[spec_RF, f, err] = dmtspec(RF,[N/sampling,W],sampling,fk,2,0.05,1);
[spec_M1F, f, err] = dmtspec(M1F,[N/sampling,W],sampling,fk,2,0.05,1);
[spec_M2F, f, err] = dmtspec(M2F,[N/sampling,W],sampling,fk,2,0.05,1);


figure;
plot(f,log10(spec_RF),'LineWidth', 2); hold on
plot(f,log10(spec_M1F),'LineWidth', 2); hold on
plot(f,log10(spec_M2F),'LineWidth', 2); hold on

legend('Rest','M1','M2')
xlabel('Frequency (Hz)')
ylabel('Log(PSD)')
xlim([0 100])
title('Filtered PSD: BandPass + Notch')
grid on



% %%%%%%%%%%%%%%%%%%
% Spectrograms
% %%%%%%%%%%%%%%%%%%

range = [10 15 20]; % trial indexes 
fk = [0 45]; % freq range
tapers = [0.3 7]; % Time resolution, Freq. resoluzion 
k = floor(2*tapers(1)*tapers(2) - 1) % number of tapers used 
dn = 0.005; % sliding step
fs = 1024; % sampling 
pad = 2;

close all;
% MEAN SPECTROGRAMS for each session 
[spec_RF,f,ti] = compute_plot_spec(RF, tapers, fs, dn, fk, pad, 0.05,'Rest avg - sess 0');
[spec_M1F,f,ti2] = compute_plot_spec(M1F, tapers, fs, dn, fk, pad, 0.05,'M1 avg - sess 0');
[spec_M2F,f,ti2]  = compute_plot_spec(M2F, tapers, fs, dn, fk, pad, 0.05,'M2 avg - sess 0');

% SINGLE TRIAL SPECTROGRAMS for a few trials 
[spec_tial_RF,f,ti] = compute_plot_spec_i(RF, tapers, fs, dn, fk, pad, 0.05,'Rest - sess 0',range);
[spec_tial_M1F,f,ti2] = compute_plot_spec_i(M1F, tapers, fs, dn, fk, pad, 0.05,'M1 - sess 0',range);
[spec_tial_M2F,f,ti2] = compute_plot_spec_i(M2F, tapers, fs, dn, fk, pad, 0.05,'M2 - sess 0',range);

 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SAVE SPECTROGRAMS TO BE USED FOR CLASSIFICATION IN PYTHON 

% Compute trail-by-trial spectrograms, for each channel separately, and save them into Data folder
[specRch_tot,f,ti] = compute_stack_save_spec_ch(RFch, tapers, fs, dn, fk, pad, 0.05, sprintf('Rspec_sess%d_ch_split_ord.mat',sess),dir_data);
[specM1ch_tot,f,ti] = compute_stack_save_spec_ch(M1ch, tapers, fs, dn, fk, pad, 0.05, sprintf('M1spec_sess%d_ch_split_ord.mat',sess),dir_data);
[specM2ch_tot,f,ti] = compute_stack_save_spec_ch(M2ch, tapers, fs, dn, fk, pad, 0.05, sprintf('M2spec_sess%d_ch_split_ord.mat',sess),dir_data);




