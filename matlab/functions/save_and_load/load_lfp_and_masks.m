function [lfp, lfp_all, mask] = load_lfp_and_masks(BRAIN_reg_rec_dir,file_name)

display(['Loading LFP data ...'])

% low speed trials (for PSD)
load(strcat(BRAIN_reg_rec_dir,sprintf('\\lfp_epoch_low_speed%s.mat',file_name))); % channel x minute x trial x lfp

% lfp_all, all trials (both low and high speed) - (for spectrograms)
load(strcat(BRAIN_reg_rec_dir,sprintf('\\lfp_epoch_all_trials%s.mat',file_name))); % channel x minute x trial x lfp

% low/high speed masks
load(strcat(BRAIN_reg_rec_dir,sprintf('\\mask_low_high_speed%s.mat',file_name))); % channel x minute x trial x lfp


end 