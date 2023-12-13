
function plot_psd_20_min(psd, main_title, dir_rec, hpc_reg, save_flag,method)

fig = figure('Position', [0, 0, 1000, 3900]);
f = psd.f;

for min = 1:size(psd.B,1) % numb of minute
    subplot(5,4,min)
    plot(f,log10(psd.B(min,:)),'LineWidth', 2, 'Color','#417C81'); hold on
    plot(f,log10(psd.L(min,:)),'LineWidth', 2, 'Color','#7986AC'); hold on
    plot(f,log10(psd.M(min,:)),'LineWidth', 2, 'Color','#2F4172'); hold on
    plot(f,log10(psd.H(min,:)),'LineWidth', 2, 'Color','#061439'); 
    
    title(sprintf('min = %d',min))
    grid on
    set(gcf, 'DefaultAxesLooseInset', [0,0,0,0]);
    
    x_ticks_max = ceil(max(f)/25) * 25; % finds the nearest multiple of 25 that's greater than or equal to max(x)
    x_ticks = 0:25:x_ticks_max;
    xticks(x_ticks);
end

ax = gca;
% Retrieve the color order matrix
defaultColors = ax.ColorOrder;

% Display the default colors
disp(defaultColors);

h_legend = legend('base', 'low','mid','high'); 
set(h_legend, 'Position', [0.92 0.85 0.07 0.03])

% Title for the entire figure (the manual way without suptitle)
fig_title = uicontrol('Style', 'text',...
    'String', main_title,...
    'Units', 'normalized',...
    'Position', [0.3 0.95 0.4 0.04],...
    'BackgroundColor', get(gcf, 'Color'),...
    'FontSize', 12,...
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