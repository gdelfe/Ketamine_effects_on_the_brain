
% Aggregate trials relative to the same minute, to compute PSD for each
% specific minute. Trials are

function [lfp_agg] = aggregate_trials_ch_range(lfp,ch_range)

lfp_agg.B = cell(size(lfp.B,2),1); % cell num of min
lfp_agg.L = cell(size(lfp.L,2),1); % cell num of minx
lfp_agg.M = cell(size(lfp.M,2),1); % cell num of min
lfp_agg.H = cell(size(lfp.H,2),1); % cell num of min

% aggregate trials for different channels together, minute by minute 
for min = 1:size(lfp.B,2) % for all the min 
    for ch = ch_range % for all the channels 
        lfp_agg.B{min} = cat(1,lfp_agg.B{min},lfp.B{ch,min}); % stack along 1st dimension all the channels for the same minute, final shape: min x trial x lfp values (1 min)
        lfp_agg.L{min} = cat(1,lfp_agg.L{min},lfp.L{ch,min}); % stack along 1st dimension all the channels for the same minute
        lfp_agg.M{min} = cat(1,lfp_agg.M{min},lfp.M{ch,min}); % stack along 1st dimension all the channels for the same minute
        lfp_agg.H{min} = cat(1,lfp_agg.H{min},lfp.H{ch,min}); % stack along 1st dimension all the channels for the same minute
    end
end


end