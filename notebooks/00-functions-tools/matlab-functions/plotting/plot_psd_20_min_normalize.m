% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot PSD for a given HPC stratum or the whole HPC depending on the chosen
% input by the user "hpc_reg". PSD is plotted for 20 min, one minute at the
% time. Both baseline, low, mid, high dose PSD are plotted together.
% PSD is normalized 
%
% Inputs: psd structure, dir_rec: path for the saving directory; hpc_reg:
% hpc region to plot; save_flag: 0/1 depending if want to not save(save);
% method: method used to appear in title; color: assigned color (don't need
% to change
%
% Output: Figures saved in given directory .../psd
%
% @ Gino Del Ferraro, Fenton Lab, Nov/Dec 2023
function plot_psd_20_min_normalize(psd, main_title, dir_rec, hpc_reg, save_flag, method, color)

fig = figure('Position', [0, 0, 700, 700]);
f = psd.f;
delta = 0.25; % shift yrange for plotting

% find min and max values across all the minutes for range in plotting
min_B = min(min(log10(psd.B(:,:)./sum(psd.B(:,:),2))));
min_L = min(min(log10(psd.L(:,:)./sum(psd.L(:,:),2))));
min_M = min(min(log10(psd.M(:,:)./sum(psd.M(:,:),2))));
min_H = min(min(log10(psd.H(:,:)./sum(psd.H(:,:),2))));

max_B = max(max(log10(psd.B(:,:)./sum(psd.B(:,:),2))));
max_L = max(max(log10(psd.L(:,:)./sum(psd.L(:,:),2))));
max_M = max(max(log10(psd.M(:,:)./sum(psd.M(:,:),2))));
max_H = max(max(log10(psd.H(:,:)./sum(psd.H(:,:),2))));


min_tot = min([min_B,min_L,min_M,min_H]);
max_tot = max([max_B,max_L,max_M,max_H]) + delta;

ha = tight_subplot(5,4,[.021 .02],[.05 .05],[.07 .2]);  % [gap_h gap_w] , margin [lower upper] [left right]
for minute = 1:size(psd.B,1) % numb of minute
    %     ax = subplot(5,4,minute);
    axes(ha(minute));
    plot(f,log10(psd.B(minute,:)/sum(psd.B(minute,:))),'LineWidth', 2, 'Color',color.B); hold on
    plot(f,log10(psd.L(minute,:)/sum(psd.L(minute,:))),'LineWidth', 2, 'Color',color.L); hold on
    plot(f,log10(psd.M(minute,:)/sum(psd.M(minute,:))),'LineWidth', 2, 'Color',color.M); hold on
    plot(f,log10(psd.H(minute,:)/sum(psd.H(minute,:))),'LineWidth', 2, 'Color',color.H);
    
    ax = gca;
    ax.LineWidth = 1;
    box off
    %title(sprintf('min = %d',minute))
    grid off;
    
    x_ticks_max = ceil(max(f)/25) * 25; % finds the nearest multiple of 25 that's greater than or equal to max(x)
    x_ticks = 0:25:x_ticks_max;
    xticks(x_ticks);
    ylim([min_tot,max_tot]);
    
    if minute == 1
        xlabel('Frequency (Hz)')
        ylabel('Norm log(PSD)')
    end
    % Remove y-axis and x-axis numbers except on the first subplot
    if minute ~= 1
        yticklabels([]);
        xticklabels([]);
    else
        xticks(x_ticks);
    end
    
    h_dummy = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
    h = legend(h_dummy,sprintf('min = %d',minute + 5));
    set(h, 'EdgeColor', 'none', 'FontSize', 10, 'Location', 'northeast');
    h.Color ='none';
    
end

h_legend = legend('base', 'low','mid','high');
set(h_legend, 'Position', [0.83 0.88 0.07 0.03])

% Title for the entire figure (the manual way without suptitle)
fig_title = uicontrol('Style', 'text',...
    'String', main_title,...
    'Units', 'normalized',...
    'Position', [0.3 0.95 0.4 0.04],...
    'BackgroundColor', get(gcf, 'Color'),...
    'FontSize', 9,...
    'FontWeight', 'bold');


if save_flag
    % Saving
    dir_out = strcat(dir_rec,['\Figures\psd\',method]);
    if ~exist(dir_out, 'dir')
        mkdir(dir_out)
    end
    saveas(fig, strcat(dir_out,sprintf('\\20_min_psd_normalized_%s.png',hpc_reg) ))
    saveas(fig, strcat(dir_out,sprintf('\\20_min_psd_normalized_%s.pdf',hpc_reg) ))
end


end