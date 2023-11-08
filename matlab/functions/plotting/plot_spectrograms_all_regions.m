% Plot spectrograms for each HPC brain subregion for a given minute of
% data, for a given epoch (base, low dose, mid dose, high dose)
%
% @ Gino Del Ferraro, NYU, Oct 2023

function plot_spectrograms_all_regions(spec_rec, epoch, mask, step_t, step_f, min, min_lab, main_title, dir_rec, save)

epoch_low = [epoch,'_low'];
epoch_high = [epoch,'_high'];
    
ts = spec_rec.ts;
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

% CA1
subplot(6,1,1)
tvimage(zscore(log10(spec_rec.(epoch).CA1(:,:,min)),1,2)); colorbar; hold on 
title(sprintf('CA1, min = %d',min_lab),'FontSize',12)
set(gca, 'XTick',valx_idx, 'XTickLabel',round(valxlbl))
set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
grid on
high_trials = high_trial_times(mask,epoch_low,epoch_high,ti,min_lab); % get high trials time 

for ii = 1:length(high_trials)
    if high_trials(ii) == 1
        plot([ii-0.5, ii+0.5], [20,20], 'k-', 'LineWidth', 1);
    end
end
 

% Ripple
subplot(6,1,2)
tvimage(zscore(log10(spec_rec.(epoch).ripple(:,:,min)),1,2)); colorbar; hold on 
title(sprintf('Ripple, min = %d',min_lab),'FontSize',12)
set(gca, 'XTick',valx_idx, 'XTickLabel',round(valxlbl))
set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
grid on

for ii = 1:length(high_trials)
    if high_trials(ii) == 1
        plot([ii-0.5, ii+0.5], [20,20], 'k-', 'LineWidth', 1);
    end
end



% Radiatum
subplot(6,1,3)
tvimage(zscore(log10(spec_rec.(epoch).rad(:,:,min)),1,2)); colorbar; hold on 
title(sprintf('Radiatum, min = %d',min_lab),'FontSize',12)
set(gca, 'XTick',valx_idx, 'XTickLabel',round(valxlbl))
set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
grid on

for ii = 1:length(high_trials)
    if high_trials(ii) == 1
        plot([ii-0.5, ii+0.5], [20,20], 'k-', 'LineWidth', 1);
    end
end


% LocMol
subplot(6,1,4)
tvimage(zscore(log10(spec_rec.(epoch).lm(:,:,min)),1,2)); colorbar; hold on 
title(sprintf('LocMol, min = %d',min_lab),'FontSize',12)
set(gca, 'XTick',valx_idx, 'XTickLabel',round(valxlbl))
set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
grid on

for ii = 1:length(high_trials)
    if high_trials(ii) == 1
        plot([ii-0.5, ii+0.5], [20,20], 'k-', 'LineWidth', 1);
    end
end


% Dentate Up
subplot(6,1,5)
tvimage(zscore(log10(spec_rec.(epoch).dup(:,:,min)),1,2)); colorbar; hold on 
title(sprintf('Dentate Up, min = %d',min_lab),'FontSize',12)
set(gca, 'XTick',valx_idx, 'XTickLabel',round(valxlbl))
set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
grid on

for ii = 1:length(high_trials)
    if high_trials(ii) == 1
        plot([ii-0.5, ii+0.5], [20,20], 'k-', 'LineWidth', 1);
    end
end

% Dentate Down
subplot(6,1,6)
tvimage(zscore(log10(spec_rec.(epoch).dd(:,:,min)),1,2)); colorbar; hold on 
title(sprintf('Dentate Down, min = %d',min_lab),'FontSize',12)
set(gca, 'XTick',valx_idx, 'XTickLabel',round(valxlbl))
set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
grid on

for ii = 1:length(high_trials)
    if high_trials(ii) == 1
        plot([ii-0.5, ii+0.5], [20,20], 'k-', 'LineWidth', 1);
    end
end


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
    dir_out = strcat(dir_rec,'\Figures\spectrograms\sub_regions\');
    if ~exist(dir_out, 'dir')
        mkdir(dir_out)
    end
    saveas(fig,strcat(dir_out,sprintf('\\spec_all_regions_epoch_%s_min_%d.png',epoch,min_lab) ) )
end

end