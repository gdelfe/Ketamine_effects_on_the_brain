
function plot_psd_all_HPC_region_all_epochs_compact(psd, dir_rec, min, save_flag, method)

fig = figure('Position', [0, 0, 800, 280]);
f = psd.so.f;
fmax = 100;

min_tot = 5.5;
max_tot = 9;
epoch = ['B','L','M','H'];
ep_title = {'Baseline', 'Low Dose', 'High Dose', 'Mid Dose'};

ha = tight_subplot(1,4,[.03 .02],[.15 .2],[.06 .07]); % [gap_h gap_w] , margin [lower upper] [left right]

for ii = 1:4
    % BASELINE
    axes(ha(ii));
    plot(f,log10(psd.so.(epoch(ii))(min,:)),'LineWidth', 2, 'Color', '#ef476f'); hold on
    plot(f,log10(psd.sp.(epoch(ii))(min,:)),'LineWidth', 2, 'Color', '#f78c6b'); hold on
    plot(f,log10(psd.rad.(epoch(ii))(min,:)),'LineWidth', 2, 'Color', '#ffd166'); hold on
    plot(f,log10(psd.lm.(epoch(ii))(min,:)),'LineWidth', 2, 'Color', '#06d6a0'); hold on
    plot(f,log10(psd.dup.(epoch(ii))(min,:)),'LineWidth', 2, 'Color', '#118ab2'); hold on
    plot(f,log10(psd.dd.(epoch(ii))(min,:)),'LineWidth', 2, 'Color', '#073b4c');
    
    ax = gca;
    ax.LineWidth = 1.2;
    ax.TickLength = [0.02, 0.025]
    box off
    %title(sprintf('min = %d',minute))
    xticks([0, 25,50,75,100])
    
    grid off;
    title(ep_title{ii},'FontSize',9, 'FontWeight','light')
    if ii == 1
        xlabel('Frequency (Hz)')
        ylabel('Norm log(PSD)')
    end
    if ii == 4
        h = legend('Oriens', 'Pyramidale','Radiatum','LocMol','Dentate Up','Dentate Down');
        set(h, 'EdgeColor', 'none', 'FontSize', 8, 'Location', 'northeast');
        h.Color ='none';
    end
    % Remove y-axis and x-axis numbers except on the first subplot
    if ii ~= 1
        yticklabels([]);
        xticklabels([]);
    end
    
    xlim([0,fmax])
    ylim([min_tot,max_tot])
    
    
    
end




% Title for the entire figure (the manual way without suptitle)
fig_title = uicontrol('Style', 'text',...
    'String', sprintf('PSD all epochs, not normalized, min = %d',min),...
    'Units', 'normalized',...
    'Position', [0.3 0.95 0.4 0.04],...
    'BackgroundColor', get(gcf, 'Color'),...
    'FontSize', 10,...
    'FontWeight', 'bold');


if save_flag
    % Saving
    dir_out = strcat(dir_rec,['\Figures\psd\',method]);
    if ~exist(dir_out, 'dir')
        mkdir(dir_out)
    end
    saveas(fig, strcat(dir_out,sprintf('\\psd_all_epochs_all_regions_min_%d_%s.png', min, method) ))
    saveas(fig, strcat(dir_out,sprintf('\\psd_all_epochs_all_regions_min_%d_%s.pdf', min, method) ))
end


end