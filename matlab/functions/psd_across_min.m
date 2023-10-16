
function [spec_B, spec_L, spec_M, spec_H, f] = spectrograms_across_min(lfp_B_all, lfp_L_all, lfp_M_all, lfp_H_all, W, fk)

display(['Computing PSD for the whole recoding  ...'])
% W = 3; % Frequency resolution for PSD 
N = 1250; % Lenght of time series
sampling = 1250; % sampling 
% fk = 100; % max frequency

spec_B = []; spec_L = []; spec_M = []; spec_H = [];
for min =1:size(lfp_B_all,1) % for all the min 
    [spec_B_min, f, err] = dmtspec(lfp_B_all{min},[N/sampling,W],sampling,fk,2,0.05,1);
    spec_B = [spec_B; spec_B_min];
    [spec_L_min, f, err] = dmtspec(lfp_L_all{min},[N/sampling,W],sampling,fk,2,0.05,1);
    spec_L = [spec_L; spec_L_min];
    [spec_M_min, f, err] = dmtspec(lfp_M_all{min},[N/sampling,W],sampling,fk,2,0.05,1);
    spec_M = [spec_M; spec_M_min];
    [spec_H_min, f, err] = dmtspec(lfp_H_all{min},[N/sampling,W],sampling,fk,2,0.05,1);
    spec_H = [spec_H; spec_H_min];
end

end 