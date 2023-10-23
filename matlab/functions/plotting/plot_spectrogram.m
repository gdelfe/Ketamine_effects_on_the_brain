% Plot mean spectrograms across channels and trials for a given session and plot them  
% 
% @ Gino Del Ferraro, NYU, June 2023

function plot_spectrogram(X, spec, f, ti, fs, step_t, step_f, str_title, epoch, norm)

ts = linspace(0,size(X,2)/fs,size(X,2));

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

fig = figure('Position', [50, 100, 1700, 200]);

subplot(2,1,1)
% no normalization 
if isempty(norm)
    tvimage(log10(spec))
    colorbar
    str_title = ['no norm - Epoch: ', epoch,' - ', str_title];
% zscore along frequency normalization 
elseif strcmp(norm,'zscore')
    tvimage(zscore(log10(spec),1,2))
    colorbar
    str_title = ['zscored - Epoch: ', epoch,' - ', str_title];
end

set(gca, 'XTick',valx_idx, 'XTickLabel',valxlbl)
set(gca, 'YTick',y_idx, 'YTickLabel',round(ylbl))
ylim([1 ,max(y_idx)])
xlabel('time (sec)')
ylabel('frequency (Hz)')
title(str_title)
grid on 

figure;
plot(tf_gamma,'LineWidth', 2); hold on
subplot(2,1,1)

tf_gamma = mean(log10(spec(:,f_gamma)),2);


tf_gamma = mean(zscore(log10(spec(:,f_gamma)),1,2),2);
plot(f,log10(spec_B(min,:)/sum(spec_B(min,:))),'LineWidth', 2); hold on

ax = gca; % Get the current axes handle
ax.Position = [0.1 0.1 0.8 0.8]; % [left, bottom, width, height]


end