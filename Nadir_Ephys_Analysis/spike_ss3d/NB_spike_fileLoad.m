% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% August 2016

% Adapted from NLfileload.m

% Create arrays for timestamps, cell numbers (groups in SS3D), and samples
[TT1_Timestamps, TT1_CellNumbers, TT1_Samples] = Nlx2MatSpike('TT1_cells.NTT', [1 0 1 0 1], 0, 1, []);

% Create arrays for event timestamps and TTL information
[EV_Timestamps, EventIDs, EV_TTLs, EventStrings] = Nlx2MatEV('Events.nev', [1 1 1 0 1], 0, 1, []);

% Convert 3D TT1 data into 2D data (only use 1 of 4 tetrode channels)
SamplesCSC1 = squeeze(TT1_Samples(:,1,:));

% Create commonly used variables
numSpikes=length(TT1_CellNumbers);
numClusters=length(unique(TT1_CellNumbers));
sampleRate = 32000; % Sample rate can also come from importing CSC information

% Create data file
spikeData = struct('TimestampsTT1', TT1_Timestamps, 'CellNumbersTT1', TT1_CellNumbers, 'SamplesTT1', TT1_Samples, 'EV_Timestamps', EV_Timestamps,  'EV_TTLs', EV_TTLs, 'SamplesCSC1', SamplesCSC1, 'numSpikes', numSpikes, 'numClusters', numClusters, 'sampleRate', sampleRate);
save spikeData spikeData;
