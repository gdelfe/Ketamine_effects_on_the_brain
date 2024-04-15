
function plot_psd_20_min_pre_post(psd_1, psd_2, main_title, dir_rec, hpc_reg, save_flag)


f = psd_1.f;

epochs = ["B","L","M","H"];
ep_names = ["baseline","low dose","mid dose","high dose"];

baseColor = zeros(4,3);

baseColor(1,:) = [0    0.4470    0.7410];
baseColor(2,:) = [0.8500    0.3250    0.0980];
baseColor(3,:) = [0.9290    0.6940    0.1250] ;
baseColor(4,:) = [0.4940    0.1840    0.5560];

psd_1_tot = []; psd_2_tot = [] ;

i = 1;
for epoch = epochs
    fig = figure('Position', [0, 0, 700, 700]);
    main_title_epoch = [main_title + " -- " + ep_names(i)];
    
    
    darkColor = darkenColor(baseColor(i,:), 0.5); % Darken by 50%
    ha = tight_subplot(5,4,[.021 .02],[.05 .05],[.07 .2]);  % [gap_h gap_w] , margin [lower upper] [left right]
    
    for minute = 1:size(psd_1.B,1) % numb of minute
        axes(ha(minute));
        
        plot(f,log10(psd_1.(epoch)(minute,:)),'LineWidth', 2, 'Color', baseColor(i,:)); hold on % first day
        plot(f,log10(psd_2.(epoch)(minute,:)),'LineWidth', 2, 'Color', darkColor); hold on % after a few days (boost)
        
        % accumulate values at each minute to compute average results
        psd_1_tot = [psd_1_tot; log10(psd_1.(epoch)(minute,:))];
        psd_2_tot = [psd_2_tot; log10(psd_2.(epoch)(minute,:))];
        
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
    
    h_legend = legend('PRE','POST');
    set(h_legend, 'Position', [0.83 0.88 0.07 0.03])
    
    
    % Title for the entire figure (the manual way without suptitle)
    fig_title = uicontrol('Style', 'text',...
        'String', main_title_epoch,...
        'Units', 'normalized',...
        'Position', [0.3 0.95 0.4 0.04],...
        'BackgroundColor', get(gcf, 'Color'),...
        'FontSize', 12,...
        'FontWeight', 'bold');
    
    psd_1_avg = mean(psd_1_tot,1);
    psd_1_sem = std(psd_1_tot)/sqrt(size(psd_1_tot,1));
    psd_2_avg = mean(psd_2_tot,1);
    psd_2_sem = std(psd_2_tot)/sqrt(size(psd_2_tot,1));
    
    fig_avg = figure('Position', [300, 300, 300, 250]);
    shadedErrorBar(f,psd_1_avg, psd_1_sem,'lineprops',{'color', baseColor(i,:)},'patchSaturation',0.4); hold on
    shadedErrorBar(f,psd_2_avg, psd_2_sem,'lineprops',{'color', darkColor},'patchSaturation',0.4);
    xlabel('Frequency (Hz)')
    ylabel('Norm log(PSD)')
    legend('PRE','POST')
    title(sprintf('Avg across minutes - %s', ep_names(i)),'FontSize',10);
    
    
    if save_flag
        % Saving
        dir_out = strcat(dir_rec,['\Figures\pre_and_post\psd\']);
        if ~exist(dir_out, 'dir')
            mkdir(dir_out)
        end
        saveas(fig, strcat(dir_out,sprintf('\\20_min_psd_pre_post_epoch_%s_reg_%s.png',epoch, hpc_reg) ))
        saveas(fig_avg, strcat(dir_out,sprintf('\\20_min_psd_AVG_pre_post_epoch_%s_reg_%s.png',epoch, hpc_reg) ))
        
        saveas(fig, strcat(dir_out,sprintf('\\20_min_psd_pre_post_epoch_%s_reg_%s.pdf',epoch, hpc_reg) ))
        saveas(fig_avg, strcat(dir_out,sprintf('\\20_min_psd_AVG_pre_post_epoch_%s_reg_%s.pdf',epoch, hpc_reg) ))
    end
    
    i = i + 1;
    
    
end


end