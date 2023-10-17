
% EMG extraction:

% 300-600Hz Coherence on special channels
% estimated on a 0.5s interval every 0.5s.

MONKEYDIR = '/vol/sas2b/Goose_Multiscale_M1_Wireless';
night = '180327';
rec = '003';
drivename = 'LM1_ECOG_3';
% night = '180322';
% rec = '001';
% drivename = 'LM1_ECOG_3';

rawecog_filename = [MONKEYDIR '/' night '/' rec '/rec' rec '.' drivename '.raw.dat'];
lfpecog_filename = [MONKEYDIR '/' night '/' rec '/rec' rec '.' drivename '.lfp.dat'];

expdefn_filename = [MONKEYDIR '/' night '/' rec '/rec' rec '.experiment.mat'];

load(expdefn_filename);
Nch = length(experiment.hardware.microdrive(1).electrodes);
Fs = experiment.hardware.acquisition.samplingrate;
LfpFs =1e3;
T = 100;
fmin = 300; fmax = 500;
W=100;
%
% % Load Raw ECoG for Coherence
% fid = fopen(rawecog_filename, 'r');
% data = fread(fid,[Nch,T*Fs],'short');
%
% % Select special channels
% specChans = [20,50];
% X = data(specChans(1),:);
% Y = data(specChans(2),:);
% % Compute coherence
% tic
% [coh_003,f,s1,s2]= tfcoh(X,Y,[.5,W],3e4,0.5,2e3,[],[],11,1,0);
% toc

recs = dayrecs(night);
nRecs = length(recs);
totSpec = [];
for iRec = 1:nRecs
    tic
    rec = recs{iRec}
    lfpecog_filename = [MONKEYDIR '/' night '/' rec '/rec' rec '.' drivename '.lfp.dat'];
    
    % Load Lfp-filtered ECoG for Spectrum
    fid = fopen(lfpecog_filename, 'r');
    lfp = fread(fid,[Nch,inf],'float');
    fclose(fid);
    lfp = lfp';
    
    % compute spectrogram for SW analysis
    numhistbins = 21;
    
    numfreqs = 100;
    freqlist = logspace(0,2.45,numfreqs);
    window = 10;
    noverlap = 9;
    window = window*LfpFs;
    noverlap = noverlap*LfpFs;
    
    [spec,f]= tfspec(lfp(:,50)',[10,1],1e3,1,300,[],[],0,1);
    nFreqlist = length(freqlist);
    myFFTspec = [];
    for iFreqlist = 1:nFreqlist
        ind = find(freqlist(iFreqlist)>f,1,'last');
        myFFTspec(iFreqlist,:) = spec(:,ind);
    end
    
    
    totSpec = [totSpec myFFTspec];
end

[ztotSpec,mu,sig] = zscore(log10(totSpec)');
totz = zscore(abs(sum(ztotSpec')));
badtimes = find(totz>1);
ztotSpec(badtimes,:) = [];

figure;
tvimage(ztotSpec(1:1000,:),'CLim',[-1 2]);

