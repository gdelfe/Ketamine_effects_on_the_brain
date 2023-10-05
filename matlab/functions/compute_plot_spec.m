% Compute mean spectrograms across channels and trials for a given session and plot them  
% 
% @ Gino Del Ferraro, NYU, June 2023

function [spec,f,ti] = compute_plot_spec(X, tapers, fs, dn, fk, pad, pval,str_title)

[spec, f, ti] = tfspec(X, tapers, fs, dn, fk, pad, pval,1);
ts = linspace(0,size(X,2)/fs,size(X,2));

figure;
tvimage(squeeze(log10(spec)))
colorbar

[x_idx, xlbl, y_idx, ylbl] = tfspec_labels(ts,ti,f,1,5);
[valx_idx, id] = unique(x_idx);
[valxlbl, id] = unique(xlbl);


set(gca, 'XTick',valx_idx, 'XTickLabel',valxlbl)
% set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
ylim([fk])
xlabel('time (sec)')
ylabel('frequency (Hz)')
title(str_title)
grid on 

end