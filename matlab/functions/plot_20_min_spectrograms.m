function plot_20_min_spectrograms(spec_rec,X,ti,f,fs,step_t,step_f)

ts = linspace(0,size(X,2)/fs,size(X,2));

[x_idx, xlbl, y_idx, ylbl] = tfspec_labels(ts,ti,f,step_t,step_f);
[valx_idx, id] = unique(x_idx);
[valxlbl, id] = unique(xlbl);


figure('Position', [0, 0, 1700, 3900]);

Nplots = 10 
for min = 1:Nplots % numb of minute
    subplot(Nplots,1,min)
    tvimage(zscore(log10(sq(spec_rec(:,:,1))),1,2))
    colorbar
    title(sprintf('min = %d',min),'FontSize',8)
    set(gca, 'XTick',valx_idx, 'XTickLabel',round(valxlbl))
    set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
    ylim([1 ,max(y_idx)])
    grid on
end
xlabel('time (sec)')

end