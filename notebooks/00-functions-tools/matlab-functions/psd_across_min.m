
function [psd] = spectrograms_across_min(lfp_agg, W, fk)

display(['Computing PSD for the whole recording  ...'])
% W = 3; % Frequency resolution for PSD 
N = 1250; % Lenght of time series
sampling = 1250; % sampling 
% fk = 100; % max frequency

psd_B = []; psd_L = []; psd_M = []; psd_H = [];
for min =1:size(lfp_agg.B,1) % for all the min 
    [psd_B_min, f, err] = dmtspec(lfp_agg.B{min},[N/sampling,W],sampling,fk,2,0.05,1);
    psd_B = [psd_B; psd_B_min];
    [psd_L_min, f, err] = dmtspec(lfp_agg.L{min},[N/sampling,W],sampling,fk,2,0.05,1);
    psd_L = [psd_L; psd_L_min];
    [psd_M_min, f, err] = dmtspec(lfp_agg.M{min},[N/sampling,W],sampling,fk,2,0.05,1);
    psd_M = [psd_M; psd_M_min];
    [psd_H_min, f, err] = dmtspec(lfp_agg.H{min},[N/sampling,W],sampling,fk,2,0.05,1);
    psd_H = [psd_H; psd_H_min];
end

psd.B = psd_B;
psd.L = psd_L;
psd.M = psd_M;
psd.H = psd_H;
psd.f = f;


end 