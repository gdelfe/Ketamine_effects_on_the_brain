
function plot_psd_20_min(spec_B, spec_L, spec_M, spec_H, f, main_title)

figure('Position', [0, 0, 1000, 3900]);

for min = 1:20
    subplot(5,4,min)
    plot(f,log10(spec_B(min,:)),'LineWidth', 2); hold on
    plot(f,log10(spec_L(min,:)),'LineWidth', 2); hold on
    plot(f,log10(spec_M(min,:)),'LineWidth', 2); hold on
    plot(f,log10(spec_H(min,:)),'LineWidth', 2); 
    title(sprintf('min = %d',min))
    grid on
    % Adjust subplot spacing for minimal space between them
    set(gcf, 'DefaultAxesLooseInset', [0,0,0,0]);
    
    x_ticks_max = ceil(max(f)/25) * 25; % finds the nearest multiple of 25 that's greater than or equal to max(x)
    x_ticks = 0:25:x_ticks_max;
    xticks(x_ticks);
    
end

h_legend = legend('base', 'low','mid','high'); 
set(h_legend, 'Position', [0.92 0.85 0.07 0.03]);

% Title for the entire figure (the manual way without suptitle)
fig_title = uicontrol('Style', 'text',...
    'String', main_title,...
    'Units', 'normalized',...
    'Position', [0.3 0.95 0.4 0.04],...
    'BackgroundColor', get(gcf, 'Color'),...
    'FontSize', 12,...
    'FontWeight', 'bold');

end 