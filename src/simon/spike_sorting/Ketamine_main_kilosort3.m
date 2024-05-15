%% you need to change most of the paths in this block

% addpath(genpath('D:\GitHub\Kilosort2')) % path to kilosort folder
% addpath('D:\GitHub\npy-matlab') % for converting to Phy
addpath(genpath('C:\Users\fentonadmin\Desktop\Kilosort30'))
addpath(genpath('C:\Users\fentonadmin\Desktop\Kilosort30\npy-matlab-master'))

%rootZ = 'G:\Spikes\Sample'; % the raw data binary file is in this folder
%input data
rootZ = 'Z:\NeuroPix\Ketamine\2022-06-10_14-50-33_M010_SAL_mPFC_HPC_0.3_0.3_0.3\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-100.2\';
%rootZ = 'C:\Users\fentonadmin\Desktop\Neuropixels_data_config\EunHye_data_in\2020-10-23_16-56-14\Record Node 106\experiment1\recording1\continuous\Neuropix-PXI-101.0\';
%rootZ = 'C:\Users\fentonadmin\Desktop\Neuropixels_data_config\Griffin_IN\';
%output data
% outputDir = 'C:\Users\fentonadmin\Desktop\Neuropixels_data_config\EunHye_data_out\2020-10-23_16-56-14\';
outputDir = 'D:\Ketamine\KSOUT\2022-06-10_14-50-33_M010_SAL_mPFC_HPC_0.3_0.3_0.3\mPFC\'
drOUT = outputDir;
if exist(drOUT,'dir') ==0; mkdir(drOUT);end
% rootH = 'H:\'; % path to temporary binary file (same size as data, should be on fast SSD)
rootH = 'C:\Users\fentonadmin\Desktop\Kilosort3\temp\';

pathToYourConfigFile = 'C:\Users\fentonlab\Desktop\Kilosort30\configFiles'; % take from Github folder and put it somewhere else (together with the master_file)
chanMapFile = 'neuropixPhase3B2_kilosortChanMap.mat';

ops.trange    = [0 Inf]; % time range to sort
ops.NchanTOT  = 384; % total number of channels in your recording

run(fullfile(pathToYourConfigFile, 'configFile384.m'))
ops.fproc   = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
ops.chanMap = fullfile(pathToYourConfigFile, chanMapFile);

%% this block runs all the steps of the algorithm
fprintf('Looking for data inside %s \n', rootZ)

% main parameter changes from Kilosort2 to v2.5
ops.sig        = 20;  % spatial smoothness constant for registration
ops.fshigh     = 300; % high-pass more aggresively
%ops.nblocks    = 5; % blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option. 
ops.nblocks    = 0;
% main parameter changes from Kilosort2.5 to v3.0
ops.Th       = [9 9];

% is there a channel map file in this folder?
fs = dir(fullfile(rootZ, 'chan*.mat'));
if ~isempty(fs)
    ops.chanMap = fullfile(rootZ, fs(1).name);
end

% find the binary file
fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
ops.fbinary = fullfile(rootZ, fs(1).name);

rez                = preprocessDataSub(ops);
rez                = datashift2(rez, 1);

[rez, st3, tF]     = extract_spikes(rez);

rez                = template_learning(rez, tF, st3);

[rez, st3, tF]     = trackAndSort(rez);

rez                = final_clustering(rez, tF, st3);

rez                = find_merges(rez, 1);

rootZ = fullfile(rootZ, 'kilosort3');
%mkdir(rootZ)
rezToPhy2(rez, drOUT);

%% %% if you want to save the results to a Matlab file.
% discard features in final rez file (too slow to save)
rez.cProj = [];
rez.cProjPC = [];

% final time sorting of spikes, for apps that use st3 directly
[~, isort]   = sortrows(rez.st3);
rez.st3      = rez.st3(isort, :);

% Ensure all GPU arrays are transferred to CPU side before saving to .mat
rez_fields = fieldnames(rez);
for i = 1:numel(rez_fields)
    field_name = rez_fields{i};
    if(isa(rez.(field_name), 'gpuArray'))
        rez.(field_name) = gather(rez.(field_name));
    end
end

fprintf('Saving final results in rez2  \n')
fname = fullfile(rootOUT, 'rez2.mat');
save(fname, 'rez', '-v7.3');