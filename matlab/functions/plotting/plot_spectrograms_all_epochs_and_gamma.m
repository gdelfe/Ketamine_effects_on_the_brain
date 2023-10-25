% Plot mean spectrograms across channels and trials for a given session and plot them
%
% @ Gino Del Ferraro, NYU, June 2023

function plot_spectrograms_all_epochs_and_gamma(spec_rec, mask, step_t, step_f, min, main_title, dir_rec, hpc_area, save)


ts = spec_rec.ts;
f = spec_rec.f;
ti = spec_rec.t ;
f_gamma = find(f>20 & f<50);


% generate labels for spectrogram
[x_idx, xlbl, y_idx, ylbl] = tfspec_labels(ts,ti,f,step_t,step_f);
[valx_idx, ~] = unique(x_idx);
[valxlbl, ~] = unique(xlbl);

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
subplot(8,1,1)
tvimage(zscore(log10(spec_rec.B(:,:,min)),1,2)); colorbar; 
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

subplot(8,1,2)
plot(mean(zscore(log10(spec_rec.B(:,f_gamma,min)),1,2),2))
grid on 

% LOW INJECTION 
subplot(4*2,1,3)
tvimage(zscore(log10(spec_rec.L(:,:,min)),1,2)); colorbar; 
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

subplot(8,1,4)
plot(mean(zscore(log10(spec_rec.L(:,f_gamma,min)),1,2),2))
grid on 


% MID INJECTION 

subplot(8,1,5)
tvimage(zscore(log10(spec_rec.M(:,:,min)),1,2)); colorbar; 
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

subplot(8,1,6)
plot(mean(zscore(log10(spec_rec.M(:,f_gamma,min)),1,2),2))
grid on 

% HIGH INJECTION 
subplot(8,1,7)
tvimage(zscore(log10(spec_rec.H(:,:,min)),1,2)); colorbar; hold on 
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

subplot(8,1,8)
plot(mean(zscore(log10(spec_rec.H(:,f_gamma,min)),1,2),2))
grid on 

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
    saveas(fig,strcat(dir_out,sprintf('\\spec_all_epochs_min_%d_and_gamma_%s.png',min, hpc_area) ) )
end


end