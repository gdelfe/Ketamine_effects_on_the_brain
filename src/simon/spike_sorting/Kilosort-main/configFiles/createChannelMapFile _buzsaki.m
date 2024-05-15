%  create a channel map file

%% linear array 1x32
Nchannels = 32;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;
xcoords   = ones(Nchannels,1);
ycoords   = [1:Nchannels]';
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)

fs = 30000; % sampling frequency
save('chanMap_1x32.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')

%% 4x8 linear
Nchannels = 32;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;

xcoords   = repmat([1 2 3 4]', 1, Nchannels/4);
xcoords   = xcoords(:);
ycoords   = repmat(1:Nchannels/4, 4, 1);
ycoords   = ycoords(:);
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)

plot(xcoords,ycoords,'k.','MarkerSize',20)

fs = 30000; % sampling frequency

save('C:\DATA\Spikes\Piroska\chanMap.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
%%

%% 4x8 buzsaki32
Nchannels = 32;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;

xcoords = [-1,-0.75,-0.5,-0.25,0, 0.25,0.5,0.75];
xcoords = xcoords + 2; %make it non-negative
xcoords = cat(2,0+xcoords,5+xcoords,10+xcoords,15+xcoords);
    
ycoords = [1 3 5 7 8 6 4 2];
ycoords = 8-ycoords+1; %
ycoords = cat(2,ycoords,ycoords,ycoords,ycoords);
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)

plot(xcoords,ycoords,'k.','MarkerSize',20)
vals = 1:32';
vals = num2cell(vals);
text(xcoords+0.1,ycoords+0.1,vals)

fs = 30000; % sampling frequency
save('chanMap_buzsaki32_v2.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
%%


%4x8 buzsaki32 v3

% Makes channel map for Buzsaki probe 4x8 ( ATLAS Neuroengineering E32B-20-S4-L10-200).
% Mapping is made for ordered channels.

chanMap = [1 2 3 4 5 6 7 8 ...
9 10 11 12 13 14 15 16 ...
17 18 19 20 21 22 23 24 ...
25 26 27 28 29 30 31 32 ];

connected = true(32, 1);

xcoords = [0 -13 13 -20 20 -27 27 -34 ...
200 187 213 180 220 173 227 166 ...
400 387 413 380 420 373 427 366 ...
600 587 613 580 620 573 627 566];

ycoords = [0 20 40 60 80 100 120 140 ...
0 20 40 60 80 100 120 140 ...
0 20 40 60 80 100 120 140 ...
0 20 40 60 80 100 120 140];

ycoords = -ycoords;

kcoords = [1 1 1 1 1 1 1 1 ...
2 2 2 2 2 2 2 2 ...
3 3 3 3 3 3 3 3 ...
4 4 4 4 4 4 4 4];

plot(xcoords,ycoords,'k.','MarkerSize',20)
vals = 1:32';
vals = num2cell(vals);
text(xcoords+0.1,ycoords+0.1,vals)

fs = 30000; % sampling frequency
save('chanMap_buzsaki32_v3.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')







%4x8 buzsaki32 v3

% Makes channel map for Buzsaki probe 4x8 ( ATLAS Neuroengineering E32B-20-S4-L10-200).
% Mapping is made for ordered channels.

chanMap = [1 2 3 4 5 6 7 8 ...
9 10 11 12 13 14 15 16 ...
17 18 19 20 21 22 23 24 ...
25 26 27 28 29 30 31 32 ];

connected = true(32, 1);

xcoords = [0 -13 13 -20 20 -27 27 -34 ...
200 187 213 180 220 173 227 166 ...
400 387 413 380 420 373 427 366 ...
600 587 613 580 620 573 627 566];
xcoords = xcoords + 35;

ycoords = [0 20 40 60 80 100 120 140 ...
0 20 40 60 80 100 120 140 ...
0 20 40 60 80 100 120 140 ...
0 20 40 60 80 100 120 140];

%ycoords = -ycoords;

kcoords = [1 1 1 1 1 1 1 1 ...
2 2 2 2 2 2 2 2 ...
3 3 3 3 3 3 3 3 ...
4 4 4 4 4 4 4 4];

plot(xcoords,ycoords,'k.','MarkerSize',20)
vals = 1:32';
vals = num2cell(vals);
text(xcoords+0.1,ycoords+0.1,vals)

fs = 30000; % sampling frequency
save('chanMap_buzsaki32_v4.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')



%4x8 buzsaki32 5

% Makes channel map for Buzsaki probe 4x8 ( ATLAS Neuroengineering E32B-20-S4-L10-200).
% Mapping is made for ordered channels.

%EXCLUDE 8 channels

%40 channels, exclude 8
Nchannels = 40;

chanMap = [1 2 3 4 5 6 7 8 ...
9 10 11 12 13 14 15 16 ...
17 18 19 20 21 22 23 24 ...
25 26 27 28 29 30 31 32 ...
33 34 35 36 37 38 39 40];

chanMap0ind = chanMap - 1;

connected = cat(1,true(32, 1),false(8,1));

xcoords = [0 -13 13 -20 20 -27 27 -34 ...
200 187 213 180 220 173 227 166 ...
400 387 413 380 420 373 427 366 ...
600 587 613 580 620 573 627 566 ...
0 0 0 0 0 0 0 0];

ycoords = [0 20 40 60 80 100 120 140 ...
0 20 40 60 80 100 120 140 ...
0 20 40 60 80 100 120 140 ...
0 20 40 60 80 100 120 140 ...
0 0 0 0 0 0 0 0];

ycoords = -ycoords;

kcoords = [1 1 1 1 1 1 1 1 ...
2 2 2 2 2 2 2 2 ...
3 3 3 3 3 3 3 3 ...
4 4 4 4 4 4 4 4 ...
0 0 0 0 0 0 0 0];

plot(xcoords,ycoords,'k.','MarkerSize',20)
vals = 1:40';
vals = num2cell(vals);
text(xcoords+0.1,ycoords+0.1,vals)

fs = 30000; % sampling frequency
save('chanMap_buzsaki32_v5.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')



% kcoords is used to forcefully restrict templates to channels in the same
% channel group. An option can be set in the master_file to allow a fraction 
% of all templates to span more channel groups, so that they can capture shared 
% noise across all channels. This option is

% ops.criterionNoiseChannels = 0.2; 

% if this number is less than 1, it will be treated as a fraction of the total number of clusters

% if this number is larger than 1, it will be treated as the "effective
% number" of channel groups at which to set the threshold. So if a template
% occupies more than this many channel groups, it will not be restricted to
% a single channel group. 