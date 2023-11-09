
function plot_psd_all_HPC_region_all_epochs_normalized(psd, dir_rec, min, save_flag, method)

fig = figure('Position', [0, 0, 1000, 3900]);
f = psd.CA1.f;
fmax = 100;
% BASELINE 
subplot(2,2,1)
plot(f,log10(psd.CA1.B(min,:)/sum(psd.CA1.B(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.ripple.B(min,:)/sum(psd.ripple.B(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.rad.B(min,:)/sum(psd.rad.B(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.lm.B(min,:)/sum(psd.lm.B(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.dup.B(min,:)/sum(psd.dup.B(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.dd.B(min,:)/sum(psd.dd.B(min,:))),'LineWidth', 2); hold on

xlabel('Frequency (Hz)')
ylabel('Log10 PSD ')
title('Baseline','FontSize',12)
xlim([0,fmax])

h_legend = legend('CA1', 'Ripple','Radiatum','LocMol','Dentate Up','Dentate Down'); 
set(h_legend)

% LOW DOSE
subplot(2,2,2)
plot(f,log10(psd.CA1.L(min,:)/sum(psd.CA1.L(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.ripple.L(min,:)/sum(psd.ripple.L(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.rad.L(min,:)/sum(psd.rad.L(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.lm.L(min,:)/sum(psd.lm.L(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.dup.L(min,:)/sum(psd.dup.L(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.dd.L(min,:)/sum(psd.dd.L(min,:))),'LineWidth', 2); hold on

xlabel('Frequency (Hz)')
ylabel('Log10 PSD ')
title('Low Dose','FontSize',12)
xlim([0,fmax])

h_legend = legend('CA1', 'Ripple','Radiatum','LocMol','Dentate Up','Dentate Down'); 
set(h_legend)


% MID DOSE
subplot(2,2,3)
plot(f,log10(psd.CA1.M(min,:)/sum(psd.CA1.M(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.ripple.M(min,:)/sum(psd.ripple.M(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.rad.M(min,:)/sum(psd.rad.M(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.lm.M(min,:)/sum(psd.lm.M(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.dup.M(min,:)/sum(psd.dup.M(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.dd.M(min,:)/sum(psd.dd.M(min,:))),'LineWidth', 2); hold on

xlabel('Frequency (Hz)')
ylabel('Log10 PSD ')
title('Mid Dose','FontSize',12)
xlim([0,fmax])

h_legend = legend('CA1', 'Ripple','Radiatum','LocMol','Dentate Up','Dentate Down'); 
set(h_legend)

% HIGH DOSE
subplot(2,2,4)
plot(f,log10(psd.CA1.H(min,:)/sum(psd.CA1.H(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.ripple.H(min,:)/sum(psd.ripple.H(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.rad.H(min,:)/sum(psd.rad.H(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.lm.H(min,:)/sum(psd.lm.H(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.dup.H(min,:)/sum(psd.dup.H(min,:))),'LineWidth', 2); hold on
plot(f,log10(psd.dd.H(min,:)/sum(psd.dd.H(min,:))),'LineWidth', 2); hold on

xlabel('Frequency (Hz)')
ylabel('Log10 PSD ')
title('High Dose','FontSize',12)
xlim([0,fmax])

h_legend = legend('CA1', 'Ripple','Radiatum','LocMol','Dentate Up','Dentate Down'); 
set(h_legend)


% Title for the entire figure (the manual way without suptitle)
fig_title = uicontrol('Style', 'text',...
    'String', sprintf('PSD all epochs, normalized, min = %d',min),...
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
    saveas(fig, strcat(dir_out,sprintf('\\psd_normalized_all_epochs_all_regions_min_%d_%s.png',min,method) ))
end


end 

