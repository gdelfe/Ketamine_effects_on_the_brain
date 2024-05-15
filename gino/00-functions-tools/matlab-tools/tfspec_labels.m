% This function creates x and y labels for a spectrogram 
% 
% INPUT: ts = time series values
%        ti = time values output of tfspec function
%        step_t = increment in the time direction for the t-labels
%        step_f = increment in the freq direction for the f-labels
%                                                         
% OUTPUT: x_idx of the label and xlbl values (same for y)
% 
% @ Gino Del Ferraro, NYU, 2023

function [x_idx, xlbl, y_idx, ylbl] = tfspec_labels(ts,ti,f,step_t,step_f)

% x labels 
x_idx = [];
val_min = round(min(ts(round(ti))),2);
val_max = round(max(ts(round(ti))),2);
val = [fliplr(0:(-step_t):val_min), 0:step_t:val_max]; % create range and ticklabels for x axis 
val = unique(val); % remove the repeated zero 
tsi = ts(round(ti));

for i = val
    [ d, ix ] = min(abs(tsi - i));
    x_idx = [x_idx, ix];
end
xlbl = round(tsi(x_idx),2);                  % New 'XTickLabel' Vector


% y labels 
y_idx = [];
val = 0:step_f:f(end);
for i=val
    [ d, ix ] = min(abs(f - i));
    y_idx = [y_idx, ix];
end
ylbl = round(f(y_idx));                  % New 'YTickLabel' Vector

end 