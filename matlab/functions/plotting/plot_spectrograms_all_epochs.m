% Plot mean spectrograms across channels and trials for a given session and plot them
%
% @ Gino Del Ferraro, NYU, June 2023

function plot_spectrograms_all_epochs(X, spec_rec, fs, step_t, step_f, min)

ts = linspace(0,size(X,2)/fs,size(X,2));
f = spec_rec.f;
ti = spec_rec.t ;

% generate labels for spectrogram
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

fig = figure('Position', [0, 0, 1700, 3900]);

subplot(4,1,1)
tvimage(zscore(log10(spec_rec.B(:,:,min)),1,2)); colorbar
title(sprintf('BASE, min = %d',min),'FontSize',12)
set(gca, 'XTick',valx_idx, 'XTickLabel',round(valxlbl))
set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
grid on

subplot(4,1,2)
tvimage(zscore(log10(spec_rec.L(:,:,min)),1,2)); colorbar
title(sprintf('LOW, min = %d',min),'FontSize',12)
set(gca, 'XTick',valx_idx, 'XTickLabel',round(valxlbl))
set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
grid on

subplot(4,1,3)
tvimage(zscore(log10(spec_rec.M(:,:,min)),1,2)); colorbar
title(sprintf('MID, min = %d',min),'FontSize',12)
set(gca, 'XTick',valx_idx, 'XTickLabel',round(valxlbl))
set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
grid on

subplot(4,1,4)
tvimage(zscore(log10(spec_rec.H(:,:,min)),1,2)); colorbar
title(sprintf('HIGH, min = %d',min),'FontSize',12)
set(gca, 'XTick',valx_idx, 'XTickLabel',round(valxlbl))
set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
grid on

xlabel('time (sec)')



end