% plot several spectrogram together, one for every min of recording

function plot_20_min_spectrograms(spec, spec_rec, step_t, step_f, minRange, epoch, main_title, dir_rec, save)

ts = spec_rec.ts;
ti = spec_rec.t;
f = spec_rec.f;

% get the labels for the spectrogram plotting 
[x_idx, xlbl, y_idx, ylbl] = tfspec_labels(ts,ti,f,step_t,step_f);
[valx_idx, ~] = unique(x_idx);
[valxlbl, ~] = unique(xlbl);

fig = figure('Position', [0, 0, 1700, 3900]);

for min = minRange % numb of minute
    subplot(minRange(end),1,min)
    tvimage(zscore(log10(sq(spec(:,:,min))),1,2))
    colorbar
    title(sprintf('%s, min = %d',epoch, min),'FontSize',8)
    set(gca, 'XTick',valx_idx, 'XTickLabel',round(valxlbl))
    set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
%     ylim([1 ,max(y_idx)])
    grid on
end

xlabel('time (sec)')


xlabel('time (sec)')

% Title for the entire figure (the manual way without suptitle)
fig_title = uicontrol('Style', 'text',...
    'String', main_title,...
    'Units', 'normalized',...
    'Position', [0.3 0.95 0.4 0.04],...
    'BackgroundColor', get(gcf, 'Color'),...
    'FontSize', 12,...
    'FontWeight', 'bold');

% Saving 
if save
    dir_out = strcat(dir_rec,'\Figures\spectrograms');
    if ~exist(dir_out, 'dir')
        mkdir(dir_out)
    end
    saveas(fig, strcat(dir_out,sprintf('\\20_min_spectrograms_%s_min_%d_%d.jpg',epoch,minRange(1),minRange(end))) ) 
end



end