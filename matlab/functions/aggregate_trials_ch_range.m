
% Aggregate trials relative to the same minute, to compute PSD for each
% specific minute. Trials are

function [lfp_agg] = aggregate_trials_ch_range(lfp_agg, lfp,reg,ch_range)

lfp_agg.B.(reg) = cell(size(lfp.B,2),1); % cell numb of min 
lfp_agg.L.(reg) = cell(size(lfp.L,2),1); % cell numb of min 
lfp_agg.M.(reg) = cell(size(lfp.M,2),1); % cell numb of min 
lfp_agg.H.(reg) = cell(size(lfp.H,2),1); % cell numb of min 

% aggregate trials for different channels together, minute by minute 
for min = 1:size(lfp.B,2) % for all the min 
    for ch = ch_range % for all the channels 
        lfp_agg.B.(reg){min} = cat(1,lfp_agg.B.(reg){min},lfp.B{ch,min}); % stack along 1st dimension all the channels for the same minute, final shape: min x trial x lfp values (1 min)
        lfp_agg.L.(reg){min} = cat(1,lfp_agg.L.(reg){min},lfp.L{ch,min}); % stack along 1st dimension all the channels for the same minute
        lfp_agg.M.(reg){min} = cat(1,lfp_agg.M.(reg){min},lfp.M{ch,min}); % stack along 1st dimension all the channels for the same minute
        lfp_agg.H.(reg){min} = cat(1,lfp_agg.H.(reg){min},lfp.H{ch,min}); % stack along 1st dimension all the channels for the same minute
    end
end


end