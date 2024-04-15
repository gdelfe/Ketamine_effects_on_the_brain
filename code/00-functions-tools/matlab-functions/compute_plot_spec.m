% Compute mean spectrograms across channels and trials for a given session and plot them  
% 
% @ Gino Del Ferraro, NYU, June 2023

function [spec,f,ti] = compute_plot_spec(X, tapers, fs, dn, fk, step_t,step_f, pad, pval,str_title)

[spec, f, ti] = tfspec(X, tapers, fs, dn, fk, pad, pval,1);
ts = linspace(0,size(X,2)/fs,size(X,2));


[x_idx, xlbl, y_idx, ylbl] = tfspec_labels(ts,ti,f,step_t,step_f);
[valx_idx, id] = unique(x_idx);
[valxlbl, id] = unique(xlbl);

% numfreqs = length(f);
% freqlist = logspace(0,2,numfreqs);
% 
% myFFTspec = [];
% for iFreqlist = 1:numfreqs
%     ind = find(freqlist(iFreqlist)>f,1,'last');
%     myFFTspec(:,iFreqlist) = spec(:,ind);
% end

figure;
tvimage(squeeze(log10(spec)))
colorbar

set(gca, 'XTick',valx_idx, 'XTickLabel',valxlbl)
set(gca, 'Yscale','log')
% set(gca, 'YTick',y_idx, 'YTickLabel',ylbl,'Yscale','log')
ylim([1 ,max(y_idx)])
xlabel('time (sec)')
ylabel('frequency (Hz)')
title(str_title)
grid on 

end