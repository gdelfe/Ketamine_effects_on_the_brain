
function [spec,f,ti] = compute_stack_save_spec_ch(X, tapers, fs, dn, fk, pad, pval,namefile,dir)

nch = 6;

spec_sess = [];
for ch = 1:nch
    % Compute spectrogram for each trials separately, FLAG = 0, channel by channel
    [spec_ch, f, ti] = tfspec(squeeze(X(ch,:,:)), tapers, fs, dn, fk, pad,pval,0);
    ts = linspace(0,size(X,2)/fs,size(X,2));
    spec_sess = cat(4,spec_sess,spec_ch); % Order: trial, time, frequency, channel 
end

spec = permute(spec_sess,[4,1,2,3]); % Order: channel, trial, time, frequency
save(strcat(dir,namefile),'spec');

end 