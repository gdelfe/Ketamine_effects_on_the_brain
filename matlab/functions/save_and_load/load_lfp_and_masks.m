function [lfp, lfp_all, mask] = load_lfp_and_masks(BRAIN_reg_rec_dir)

display(['Loading LFP data ...'])

% low speed trials (for PSD)
load(strcat(BRAIN_reg_rec_dir,'\lfp_epoch_low_speed_CAR.mat')); % channel x minute x trial x lfp

% lfp_all, all trials (both low and high speed) - (for spectrograms)
load(strcat(BRAIN_reg_rec_dir,'\lfp_epoch_all_trials_CAR.mat')); % channel x minute x trial x lfp

% low/high speed masks
load(strcat(BRAIN_reg_rec_dir,'\mask_low_high_speed_CAR.mat')); % channel x minute x trial x lfp


end 