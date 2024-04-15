
% Aggregate trials relative to the same minute, to compute PSD 

function [lfp_agg] = aggregate_trials(lfp)

lfp_agg.B = cell(size(lfp.B,2),1); % cell min x
lfp_agg.L = cell(size(lfp.L,2),1); % cell min x
lfp_agg.M = cell(size(lfp.M,2),1); % cell min x
lfp_agg.H = cell(size(lfp.H,2),1); % cell min x

% aggregate trials for different channels together, minute by minute 
for min = 1:size(lfp.B,2) % for all the min 
    for ch = 1:size(lfp.B,1) % for all the channels 
        lfp_agg.B{min} = cat(1,lfp_agg.B{min},lfp.B{ch,min}); % stack along 1st dimension all the channels for the same minute, final shape: min x trial x lfp values (1 min)
        lfp_agg.L{min} = cat(1,lfp_agg.L{min},lfp.L{ch,min}); % stack along 1st dimension all the channels for the same minute
        lfp_agg.M{min} = cat(1,lfp_agg.M{min},lfp.M{ch,min}); % stack along 1st dimension all the channels for the same minute
        lfp_agg.H{min} = cat(1,lfp_agg.H{min},lfp.H{ch,min}); % stack along 1st dimension all the channels for the same minute
    end
end


end