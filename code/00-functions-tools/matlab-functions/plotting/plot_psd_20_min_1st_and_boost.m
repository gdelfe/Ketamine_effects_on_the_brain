
function plot_psd_20_min_1st_and_boost(psd_1, psd_2, main_title, dir_rec, hpc_reg, save_flag, method)


f = psd_1.f;
epochs = ["B","L","M","H"];
ep_names = ["baseline","low dose","mid dose","high dose"];

baseColor = zeros(4,3);

baseColor(1,:) = [0    0.4470    0.7410];
baseColor(2,:) = [0.8500    0.3250    0.0980];
baseColor(3,:) = [0.9290    0.6940    0.1250] ;
baseColor(4,:) = [0.4940    0.1840    0.5560];


i = 1;
for epoch = epochs
    fig = figure('Position', [0, 0, 1000, 3900]);
    main_title_epoch = [main_title + " -- " + ep_names(i)]; 
    
    
    darkColor = darkenColor(baseColor(i,:), 0.5); % Darken by 50%
    
    for min = 1:size(psd_1.B,1) % numb of minute
        subplot(5,4,min)
        
        plot(f,log10(psd_1.(epoch)(min,:)),'LineWidth', 2, 'Color', baseColor(i,:)); hold on % first day
        plot(f,log10(psd_2.(epoch)(min,:)),'LineWidth', 2, 'Color', darkColor); hold on % after a few days (boost)
        
        title(sprintf('min = %d',min))
        grid on
        set(gcf, 'DefaultAxesLooseInset', [0,0,0,0]);
        
        x_ticks_max = ceil(max(f)/25) * 25; % finds the nearest multiple of 25 that's greater than or equal to max(x)
        x_ticks = 0:25:x_ticks_max;
        xticks(x_ticks);
    end
    
    i = i +1;
    
    h_legend = legend('first day', 'boost');
    set(h_legend, 'Position', [0.92 0.85 0.07 0.03])
    
    % Title for the entire figure (the manual way without suptitle)
    fig_title = uicontrol('Style', 'text',...
        'String', main_title_epoch,...
        'Units', 'normalized',...
        'Position', [0.3 0.95 0.4 0.04],...
        'BackgroundColor', get(gcf, 'Color'),...
        'FontSize', 12,...
        'FontWeight', 'bold');
    
    
    if save_flag
        % Saving
        dir_out = strcat(dir_rec,['\Figures\1st_and_boost\psd\',method]);
        if ~exist(dir_out, 'dir')
            mkdir(dir_out)
        end
        saveas(fig, strcat(dir_out,sprintf('\\20_min_psd_1st_and_boost_epoch_%s_reg_%s.png',epoch, hpc_reg) ))
    end
    
end


end