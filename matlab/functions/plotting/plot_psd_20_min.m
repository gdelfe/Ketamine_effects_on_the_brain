
function plot_psd_20_min(psd, main_title, dir_rec, hpc_reg, save_flag,method)

fig = figure('Position', [0, 0, 700, 700]);
f = psd.f;

delta = 0.25; % shift yrange for plotting

% find min and max values across all the minutes for range in plotting
min_B = min(min(log10(psd.B(:,:))));
min_L = min(min(log10(psd.L(:,:))));
min_M = min(min(log10(psd.M(:,:))));
min_H = min(min(log10(psd.H(:,:))));

max_B = max(max(log10(psd.B(:,:))));
max_L = max(max(log10(psd.L(:,:))));
max_M = max(max(log10(psd.M(:,:))));
max_H = max(max(log10(psd.H(:,:))));


min_tot = min([min_B,min_L,min_M,min_H]);
max_tot = max([max_B,max_L,max_M,max_H]) + delta;

ha = tight_subplot(5,4,[.021 .02],[.05 .05],[.07 .2]); 
for minute = 1:size(psd.B,1) % numb of minute
%     subplot(5,4,minute)
    axes(ha(minute));
    plot(f,log10(psd.B(minute,:)),'LineWidth', 2, 'Color','#87F07D'); hold on
    plot(f,log10(psd.L(minute,:)),'LineWidth', 2, 'Color','#6FE3F0'); hold on
    plot(f,log10(psd.M(minute,:)),'LineWidth', 2, 'Color','#589DE5'); hold on
    plot(f,log10(psd.H(minute,:)),'LineWidth', 2, 'Color','#2C45A9'); 
    
    ax = gca;
    ax.LineWidth = 1;
    box off
    %title(sprintf('min = %d',minute))
    grid off;
        
    x_ticks_max = ceil(max(f)/25) * 25; % finds the nearest multiple of 25 that's greater than or equal to max(x)
    x_ticks = 0:25:x_ticks_max;
    xticks(x_ticks);
    
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

ax = gca;
% Retrieve the color order matrix
defaultColors = ax.ColorOrder;

% Display the default colors
disp(defaultColors);

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
    saveas(fig, strcat(dir_out,sprintf('\\20_min_psd_%s.png',hpc_reg) ))
end


end 