% Plot mean spectrograms across channels and trials for a given session and plot them
%
% @ Gino Del Ferraro, NYU, June 2023

function plot_spectrograms_all_epochs(spec_rec, mask, step_t, step_f, min, main_title, dir_rec, save)


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

% BASELINE 
subplot(4,1,1)
tvimage(zscore(log10(spec_rec.B.HPC(:,:,min)),1,2)); colorbar; hold on 
title(sprintf('BASE, min = %d',min),'FontSize',12)
set(gca, 'XTick',valx_idx, 'XTickLabel',round(valxlbl))
set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
grid on
high_trials = high_trial_times(mask,'B_low','B_high',ti,min); % get high trials time 
A = spec_rec.B(:,:,min);
for ii = 1:length(high_trials)
    if high_trials(ii) == 1
        plot([ii-0.5, ii+0.5], [20,20], 'k-', 'LineWidth', 1);
    end
end
 

% LOW INJECTION 
subplot(4,1,2)
tvimage(zscore(log10(spec_rec.L.HPC(:,:,min)),1,2)); colorbar; hold on 
title(sprintf('LOW, min = %d',min),'FontSize',12)
set(gca, 'XTick',valx_idx, 'XTickLabel',round(valxlbl))
set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
grid on
high_trials = high_trial_times(mask,'L_low','L_high',ti,min); % get high trials time 
A = spec_rec.L(:,:,min);
for ii = 1:length(high_trials)
    if high_trials(ii) == 1
        plot([ii-0.5, ii+0.5], [20,20], 'k-', 'LineWidth', 1);
    end
end



% MID INJECTION 

subplot(4,1,3)
tvimage(zscore(log10(spec_rec.M.HPC(:,:,min)),1,2)); colorbar; hold on 
title(sprintf('MID, min = %d',min),'FontSize',12)
set(gca, 'XTick',valx_idx, 'XTickLabel',round(valxlbl))
set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
grid on
high_trials = high_trial_times(mask,'M_low','M_high',ti,min); % get high trials time 
A = spec_rec.M(:,:,min);
for ii = 1:length(high_trials)
    if high_trials(ii) == 1
        plot([ii-0.5, ii+0.5], [20,20], 'k-', 'LineWidth', 1);
    end
end


% HIGH INJECTION 
subplot(4,1,4)
tvimage(zscore(log10(spec_rec.H.HPC(:,:,min)),1,2)); colorbar; hold on 
title(sprintf('HIGH, min = %d',min),'FontSize',12)
set(gca, 'XTick',valx_idx, 'XTickLabel',round(valxlbl))
set(gca, 'YTick',y_idx, 'YTickLabel',ylbl)
grid on
high_trials = high_trial_times(mask,'H_low','H_high',ti,min); % get high trials time 
A = spec_rec.B(:,:,min);
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
    dir_out = strcat(dir_rec,'\Figures\spectrograms');
    if ~exist(dir_out, 'dir')
        mkdir(dir_out)
    end
    saveas(fig,strcat(dir_out,sprintf('\\spec_all_epochs_HPC_min_%d_%s.png',min) ) )
end

end