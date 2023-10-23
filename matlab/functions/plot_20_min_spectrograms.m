% plot several spectrogram together, one for every min of recording

function plot_20_min_spectrograms(spec, spec_rec, X, fs, step_t, step_f, minRange)

ts = linspace(0,size(X,2)/fs,size(X,2));
ti = spec_rec.t;
f = spec_rec.f;

% get the labels for the spectrogram plotting 
[x_idx, xlbl, y_idx, ylbl] = tfspec_labels(ts,ti,f,step_t,step_f);
[valx_idx, ~] = unique(x_idx);
[valxlbl, ~] = unique(xlbl);

figure('Position', [0, 0, 1700, 3900]);

for min = minRange % numb of minute
    subplot(minRange(end),1,min)
    tvimage(zscore(log10(sq(spec(:,:,min))),1,2))
    colorbar
    title(sprintf('min = %d',min),'FontSize',8)
    set(gca, 'XTick',valx_idx, 'XTickLabel',round(valxlbl))
    set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
%     ylim([1 ,max(y_idx)])
    grid on
end
xlabel('time (sec)')

end