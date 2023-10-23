function load_lfp_and_masks(BRAIN_reg_rec_dir)

display(['Loading LFP data ...'])
%
load(strcat(BRAIN_reg_rec_dir,'\lfp_epoch_low_speed.mat')); % channel x minute x trial x lfp
lfp.B = lfp_B
load(strcat(BRAIN_reg_rec_dir,'\lfp_L_epoch_low_speed.mat')); % channel x minute x trial x lfp
lfp.L = lfp_L
load(strcat(BRAIN_reg_rec_dir,'\lfp_M_epoch_low_speed.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\lfp_H_epoch_low_speed.mat')); % channel x minute x trial x lfp

% lfp_X_all type of variables 
load(strcat(BRAIN_reg_rec_dir,'\lfp_B_epoch_all_trials.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\lfp_L_epoch_all_trials.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\lfp_M_epoch_all_trials.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\lfp_H_epoch_all_trials.mat')); % channel x minute x trial x lfp

% low speed masks
load(strcat(BRAIN_reg_rec_dir,'\mask_low_B.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\mask_low_L.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\mask_low_M.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\mask_low_H.mat')); % channel x minute x trial x lfp

% high speed masks
load(strcat(BRAIN_reg_rec_dir,'\mask_high_B.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\mask_high_L.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\mask_high_M.mat')); % channel x minute x trial x lfp
load(strcat(BRAIN_reg_rec_dir,'\mask_high_H.mat')); % channel x minute x trial x lfp


end 