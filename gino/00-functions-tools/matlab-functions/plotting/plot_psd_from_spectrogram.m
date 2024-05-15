
% Plot PSD from spectrogram by averaging the time-dimension of the
% spectrogram 

function plot_psd_from_spectrogram(spec_rec, dir_rec, min, min_lab, save)

fig = figure('Position', [0, 0, 1700, 3900]);
f = spec_rec.f;
fmax = 100;

% BASELINE 
epoch = 'B';
psd_CA1 = mean(log10(spec_rec.(epoch).CA1(:,:,min)),1);
psd_ripple = mean(log10(spec_rec.(epoch).ripple(:,:,min)),1);
psd_rad = mean(log10(spec_rec.(epoch).rad(:,:,min)),1);
psd_lm = mean(log10(spec_rec.(epoch).lm(:,:,min)),1);
psd_dup = mean(log10(spec_rec.(epoch).dup(:,:,min)),1);
psd_dd = mean(log10(spec_rec.(epoch).dd(:,:,min)),1);

subplot(2,2,1)
plot(f,psd_CA1,'LineWidth', 2); hold on 
plot(f,psd_ripple,'LineWidth', 2); hold on 
plot(f,psd_rad,'LineWidth', 2); hold on 
plot(f,psd_lm,'LineWidth', 2); hold on 
plot(f,psd_dup,'LineWidth', 2); hold on 
plot(f,psd_dd,'LineWidth', 2); 
xlabel('Frequency (Hz)')
ylabel('Log10 PSD ')
title('Baseline','FontSize',12)
xlim([0,fmax])

h_legend = legend('CA1', 'Ripple','Radiatum','LocMol','Dentate Up','Dentate Down'); 
set(h_legend)

% LOW DOSE  
epoch = 'L';
psd_CA1 = mean(log10(spec_rec.(epoch).CA1(:,:,min)),1);
psd_ripple = mean(log10(spec_rec.(epoch).ripple(:,:,min)),1);
psd_rad = mean(log10(spec_rec.(epoch).rad(:,:,min)),1);
psd_lm = mean(log10(spec_rec.(epoch).lm(:,:,min)),1);
psd_dup = mean(log10(spec_rec.(epoch).dup(:,:,min)),1);
psd_dd = mean(log10(spec_rec.(epoch).dd(:,:,min)),1);

subplot(2,2,2)
plot(f,psd_CA1,'LineWidth', 2); hold on 
plot(f,psd_ripple,'LineWidth', 2); hold on 
plot(f,psd_rad,'LineWidth', 2); hold on 
plot(f,psd_lm,'LineWidth', 2); hold on 
plot(f,psd_dup,'LineWidth', 2); hold on 
plot(f,psd_dd,'LineWidth', 2); 
title('Low Dose','FontSize',12)
xlabel('Frequency (Hz)')
ylabel('Zscored Log PSD ')
xlim([0,fmax])

h_legend = legend('CA1', 'Ripple','Radiatum','LocMol','Dentate Up','Dentate Down'); 
set(h_legend)

% MID DOSE  
epoch = 'M';
psd_CA1 = mean(log10(spec_rec.(epoch).CA1(:,:,min)),1);
psd_ripple = mean(log10(spec_rec.(epoch).ripple(:,:,min)),1);
psd_rad = mean(log10(spec_rec.(epoch).rad(:,:,min)),1);
psd_lm = mean(log10(spec_rec.(epoch).lm(:,:,min)),1);
psd_dup = mean(log10(spec_rec.(epoch).dup(:,:,min)),1);
psd_dd = mean(log10(spec_rec.(epoch).dd(:,:,min)),1);

subplot(2,2,3)
plot(f,psd_CA1,'LineWidth', 2); hold on 
plot(f,psd_ripple,'LineWidth', 2); hold on 
plot(f,psd_rad,'LineWidth', 2); hold on 
plot(f,psd_lm,'LineWidth', 2); hold on 
plot(f,psd_dup,'LineWidth', 2); hold on 
plot(f,psd_dd,'LineWidth', 2);
title('Mid Dose','FontSize',12)
xlabel('Frequency (Hz)')
ylabel('Zscored Log PSD ')
xlim([0,fmax])

h_legend = legend('CA1', 'Ripple','Radiatum','LocMol','Dentate Up','Dentate Down'); 
set(h_legend)

% HIGH DOSE  
epoch = 'M';
psd_CA1 = mean(log10(spec_rec.(epoch).CA1(:,:,min)),1);
psd_ripple = mean(log10(spec_rec.(epoch).ripple(:,:,min)),1);
psd_rad = mean(log10(spec_rec.(epoch).rad(:,:,min)),1);
psd_lm = mean(log10(spec_rec.(epoch).lm(:,:,min)),1);
psd_dup = mean(log10(spec_rec.(epoch).dup(:,:,min)),1);
psd_dd = mean(log10(spec_rec.(epoch).dd(:,:,min)),1);

subplot(2,2,4)
plot(f,psd_CA1,'LineWidth', 2); hold on 
plot(f,psd_ripple,'LineWidth', 2); hold on 
plot(f,psd_rad,'LineWidth', 2); hold on 
plot(f,psd_lm,'LineWidth', 2); hold on 
plot(f,psd_dup,'LineWidth', 2); hold on 
plot(f,psd_dd,'LineWidth', 2);
xlim([0,fmax])
title('High Dose','FontSize',12)
xlabel('Frequency (Hz)')
ylabel('Zscored Log PSD ')

h_legend = legend('CA1', 'Ripple','Radiatum','LocMol','Dentate Up','Dentate Down'); 
set(h_legend)


% Title for the entire figure (the manual way without suptitle)
fig_title = uicontrol('Style', 'text',...
    'String', sprintf('Log10 PSD from zscored-spectrogram, min = %d',min_lab),...
    'Units', 'normalized',...
    'Position', [0.3 0.95 0.4 0.04],...
    'BackgroundColor', get(gcf, 'Color'),...
    'FontSize', 12,...
    'FontWeight', 'bold');


% Saving
if save
    dir_out = strcat(dir_rec,'\Figures\psd\sub_regions\');
    if ~exist(dir_out, 'dir')
        mkdir(dir_out)
    end
    saveas(fig,strcat(dir_out,sprintf('\\psd_from_Spec_all_epochs_all_areas_min_%d.png',min_lab) ) )
end



end 