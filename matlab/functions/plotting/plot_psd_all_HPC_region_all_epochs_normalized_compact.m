
function plot_psd_all_HPC_region_all_epochs_normalized_compact(psd, dir_rec, minute, save_flag, method)

fig = figure('Position', [0, 0, 900, 280]);
f = psd.so.f;
fmax = 100;

min_tot = -5;
max_tot = -1;
epoch = ['B','L','M','H'];
ep_title = {'Baseline', 'Low Dose', 'High Dose', 'Mid Dose'};

ha = tight_subplot(1,4,[.0 .02],[.15 .13],[.06 .07]); % [gap_h gap_w] , margin [lower upper] [left right]

for ii = 1:4
    % BASELINE
    axes(ha(ii));
      
    plot(f,log10(psd.so.(epoch(ii))(minute,:)/sum(psd.so.(epoch(ii))(minute,:))),'LineWidth', 2); hold on
    plot(f,log10(psd.sp.(epoch(ii))(minute,:)/sum(psd.sp.(epoch(ii))(minute,:))),'LineWidth', 2); hold on
    plot(f,log10(psd.rad.(epoch(ii))(minute,:)/sum(psd.rad.(epoch(ii))(minute,:))),'LineWidth', 2); hold on
    plot(f,log10(psd.lm.(epoch(ii))(minute,:)/sum(psd.lm.(epoch(ii))(minute,:))),'LineWidth', 2); hold on
    plot(f,log10(psd.dup.(epoch(ii))(minute,:)/sum(psd.dup.(epoch(ii))(minute,:))),'LineWidth', 2); hold on
    plot(f,log10(psd.dd.(epoch(ii))(minute,:)/sum(psd.dd.(epoch(ii))(minute,:))),'LineWidth', 2); hold on
    
    ax = gca;
    ax.LineWidth = 1.2;
    ax.TickLength = [0.02, 0.025];
    box off
    %title(sprintf('min = %d',minute))
    xticks([0, 25,50,75,100])
    
    grid off;
    title(ep_title{ii},'FontSize',9, 'FontWeight','light')
    if ii == 1
        h = legend('Oriens', 'Pyramidale','Radiatum','LocMol','Dentate Up','Dentate Down');
        set(h, 'EdgeColor', 'none', 'FontSize', 8, 'Location', 'northeast');
        h.Color ='none';
        xlabel('Frequency (Hz)')
        ylabel('Norm log(PSD)')
    end
    % Remove y-axis and x-axis numbers except on the first subplot
    if ii ~= 1
        yticklabels([]);
        xticklabels([]);
    end
    
%     xlim([0,fmax])
%     ylim([min_tot,max_tot])
    
    
    
end




% Title for the entire figure (the manual way without suptitle)
fig_title = uicontrol('Style', 'text',...
    'String', sprintf('PSD all epochs, normalized, min = %d',minute),...
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
    saveas(fig, strcat(dir_out,sprintf('\\psd_normalized_all_epochs_all_regions_min_%d_%s.png', minute, method) ))
end


end