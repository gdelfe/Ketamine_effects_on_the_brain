
%plot EEG from neuropixels
clear all; close all; clc;

%update path to your EEGLAB
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\');
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\sigprocfunc\'); 
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\guifunc\');
addpath('C:\Users\fentonlab\Desktop\Gilgamesh\Codes\eeglab\functions\adminfunc\');
addpath('C:\Users\fentonlab\Desktop\Gino\chronux_2_11\');
addpath('C:\Users\fentonlab\Desktop\Gino\Gino_codes\');

main_dir = 'C:\Users\fentonlab\Desktop\Gino\LFPs\HPC\2022-08-01_04-30-00_M015_RSK_mPFC_HPC_3_10_30mpk\';

load(strcat(main_dir,'lfp_B_epoch_low_speed.mat')) % channel x minute x trial x lfp
load(strcat(main_dir,'lfp_L_epoch_low_speed.mat')) % channel x minute x trial x lfp
load(strcat(main_dir,'lfp_M_epoch_low_speed.mat')) % channel x minute x trial x lfp
load(strcat(main_dir,'lfp_H_epoch_low_speed.mat')) % channel x minute x trial x lfp

figure;
plot(lfp_B{1,1}(1,:)); hold on 
plot(lfp_L{1,1}(1,:)); hold on 
plot(lfp_M{1,1}(1,:)); hold on 
plot(lfp_H{1,1}(1,:))


[lfp_B_all, lfp_L_all, lfp_M_all, lfp_H_all] = aggregate_trials(lfp_B,lfp_L,lfp_M,lfp_H)

lfp_B_all = cell(size(lfp_B,2),1); % cell min x
lfp_L_all = cell(size(lfp_L,2),1); % cell min x
lfp_M_all = cell(size(lfp_M,2),1); % cell min x
lfp_H_all = cell(size(lfp_H,2),1); % cell min x

% aggregate trials for different channels together, minute by minute 
for min = 1:size(lfp_B,2) % for all the min 
    for ch = 1:size(lfp_B,1) % for all the channels 
        lfp_B_all{min} = cat(1,lfp_B_all{min},lfp_B{ch,min}); % stack along 1st dimension all the channels for the same minute, final shape: min x trial x lfp values (1 min)
        lfp_L_all{min} = cat(1,lfp_L_all{min},lfp_L{ch,min}); % stack along 1st dimension all the channels for the same minute
        lfp_M_all{min} = cat(1,lfp_M_all{min},lfp_M{ch,min}); % stack along 1st dimension all the channels for the same minute
        lfp_H_all{min} = cat(1,lfp_H_all{min},lfp_H{ch,min}); % stack along 1st dimension all the channels for the same minute
    end
end 


W = 3; % Frequency resolution for PSD 
N = 1250; % Lenght of time series
sampling = 1250; % sampling 
fk = 100; % max frequency

spec_B = []; spec_L = []; spec_M = []; spec_H = [];
for min =1:20
    [spec_B_min, f, err] = dmtspec(lfp_B_all{min},[N/sampling,W],sampling,fk,2,0.05,1);
    spec_B = [spec_B; spec_B_min];
    [spec_L_min, f, err] = dmtspec(lfp_L_all{min},[N/sampling,W],sampling,fk,2,0.05,1);
    spec_L = [spec_L; spec_L_min];
    [spec_M_min, f, err] = dmtspec(lfp_M_all{min},[N/sampling,W],sampling,fk,2,0.05,1);
    spec_M = [spec_M; spec_M_min];
    [spec_H_min, f, err] = dmtspec(lfp_H_all{min},[N/sampling,W],sampling,fk,2,0.05,1);
    spec_H = [spec_H; spec_H_min];
end


figure;
plot(f,log10(spec_B(min,:)),'LineWidth', 2); hold on
plot(f,log10(spec_L(min,:)),'LineWidth', 2); hold on
plot(f,log10(spec_M(min,:)),'LineWidth', 2); hold on
plot(f,log10(spec_H(min,:)),'LineWidth', 2); hold on


legend('base','low','mid','high')
xlabel('Frequency (Hz)')
ylabel('Log(PSD)')
xlim([0 100])
title(sprintf('PSD for min = %d',min))
grid on

for min = 1:20
    subplot(5,4,min)
    plot(f,log10(spec_B(min,:)/sum(spec_B(min,:))),'LineWidth', 2); hold on
    plot(f,log10(spec_L(min,:)/sum(spec_L(min,:))),'LineWidth', 2); hold on
    plot(f,log10(spec_M(min,:)/sum(spec_M(min,:))),'LineWidth', 2); hold on
    plot(f,log10(spec_H(min,:)/sum(spec_H(min,:))),'LineWidth', 2); 
    title(sprintf('PSD for min = %d',min))
    grid on
    % Adjust subplot spacing for minimal space between them
    set(gcf, 'DefaultAxesLooseInset', [0,0,0,0]);
%     tight_layout();
end
% Title for the entire figure (the manual way without suptitle)
fig_title = uicontrol('Style', 'text',...
    'String', 'PSD - Stationary - RS Ketamine',...
    'Units', 'normalized',...
    'Position', [0.3 0.95 0.4 0.04],...
    'BackgroundColor', get(gcf, 'Color'),...
    'FontSize', 15,...
    'FontWeight', 'bold');




for min = 1:20
    subplot(5,4,min)
    plot(f,log10(spec_B(min,:)/sum(spec_B(min,:))),'LineWidth', 2); hold on
    plot(f,log10(spec_L(min,:)/sum(spec_L(min,:))),'LineWidth', 2); hold on
    plot(f,log10(spec_M(min,:)/sum(spec_M(min,:))),'LineWidth', 2); hold on
    plot(f,log10(spec_H(min,:)/sum(spec_H(min,:))),'LineWidth', 2); 
    title(sprintf('PSD for min = %d',min))
    grid on
    % Adjust subplot spacing for minimal space between them
    set(gcf, 'DefaultAxesLooseInset', [0,0,0,0]);
%     tight_layout();
end
% Title for the entire figure (the manual way without suptitle)
fig_title = uicontrol('Style', 'text',...
    'String', 'PDS normalized - Stationary - RS Ketamine',...
    'Units', 'normalized',...
    'Position', [0.3 0.95 0.4 0.04],...
    'BackgroundColor', get(gcf, 'Color'),...
    'FontSize', 15,...
    'FontWeight', 'bold');



