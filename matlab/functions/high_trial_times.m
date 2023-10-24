% Get time of high trials at the LFP sampling resolution 

function high_trials = high_trial_times(mask,epoch_low,epoch_high,ti,min)

L = any(sq(mask.(epoch_low)(min,:,:)),2); % low mask + artifacts
H = any(sq(mask.(epoch_high)(min,:,:)),2); % high mask + artifacts
H_trial = (H == 1) & (L == 0); % high mask 
M = round(length(ti)/60); % multiplication factor to have a mask at fs resolution 
high_trials = repmat(H_trial, [1 M]);
high_trials = reshape(high_trials', [1 numel(high_trials)]); % Hight trial mask 

end 