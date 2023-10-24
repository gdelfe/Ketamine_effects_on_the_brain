% Compute single trials spectrograms and plot them 
% 
% @ Gino Del Ferraro, NYU, June 2023


function [spec,f,ti] = compute_plot_spec_i(X, tapers, fs, dn, fk, pad, pval,str_title,range)

% Compute spectrogram for each trials separately, FLAG = 0
[spec, f, ti] = tfspec(X, tapers, fs, dn, fk, pad, 0.05,0);
ts = linspace(0,size(X,2)/fs,size(X,2));

for i= range
    
    figure;
    tvimage(squeeze(log10(spec(i,:,:))));
    colorbar
    
    [x_idx, xlbl, y_idx, ylbl] = tfspec_labels(ts,ti,f,1,5);
    [valx_idx, id] = unique(x_idx);
    [valxlbl, id] = unique(xlbl);
    
    str = append(str_title,',  trial = ',num2str([i]));
    
    set(gca, 'XTick',valx_idx, 'XTickLabel',valxlbl);
%     set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
    ylim([fk])
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
    title(str);
    grid on 
end

end