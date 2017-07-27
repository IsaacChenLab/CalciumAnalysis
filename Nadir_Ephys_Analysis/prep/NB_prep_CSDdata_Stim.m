% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% Last edited: July 5, 2016

% Converts CSC_Wave data into a 3D CSD_DataBlock
% Then, averages CSD_DataBlock to a 2D CSD_Data table
% This 2D table can be plotted as averaged CSC data

% Total period of time to be plotted
ms = 550;

% Number of miliseconds to include before the stimulation
offsetms = -50;

% Converts offsetms into a microsecond timestamp format (y-axis)
offset = (offsetms/1000)*sampleRate;

% Converts plot length into a microsecond timestamp format (x-axis)
lengthWave = (ms/1000)*sampleRate;

% Number of channels
channels = 32;

% Find all TTL events when stim is turned on
TTLon = find (EV_TTLs==1);

% Finds all timestamps when stim is turned on based on TTL events
stimOnTimestamps = EV_Timestamps(TTLon);

% Number of stims calculated by number of times stim was turned on
numStims = length(stimOnTimestamps);

% Interpolates the index of timeVector to match the timestamp of stimOnTimestamp
% interp_index is used to match stimOnTimestamp with CSC_Wave data
% This allows us to correlate stim events with CSC data
interp_index = interp1(timeVector,1:length(timeVector),stimOnTimestamps,'nearest');

% Creates CSD_Datablock with zeros based on # channels and length of plot
CSD_DataBlock = zeros (channels, lengthWave);
rateHistogramArray = zeros(1,numStims);

% Adds CSC data to CSD block for each stimulation
for thisStim = 1:numStims; 
    
    % Calculates timestamp coordinate for beginning of graph (start of stim - offset)
    startplot = interp_index(thisStim)+ offset;
    % Calculates timestamp coordinate for end of graph (why -1 at the end?)
    endplot = startplot + lengthWave - 1;
   
    % Adds each stim plot to CSD_DataBlock
    CSD_DataBlock (1,:,thisStim) = CSC1_Wave(startplot:endplot);
    CSD_DataBlock (2,:,thisStim) = CSC2_Wave(startplot:endplot);
    CSD_DataBlock (3,:,thisStim) = CSC3_Wave(startplot:endplot);
    CSD_DataBlock (4,:,thisStim) = CSC4_Wave(startplot:endplot);
    CSD_DataBlock (5,:,thisStim) = CSC5_Wave(startplot:endplot);
    CSD_DataBlock (6,:,thisStim) = CSC6_Wave(startplot:endplot);
    CSD_DataBlock (7,:,thisStim) = CSC7_Wave(startplot:endplot);
    CSD_DataBlock (8,:,thisStim) = CSC8_Wave(startplot:endplot);
    CSD_DataBlock (9,:,thisStim) = CSC9_Wave(startplot:endplot);
    CSD_DataBlock (10,:,thisStim) = CSC10_Wave(startplot:endplot);
    CSD_DataBlock (11,:,thisStim) = CSC11_Wave(startplot:endplot);
    CSD_DataBlock (12,:,thisStim) = CSC12_Wave(startplot:endplot);
    CSD_DataBlock (13,:,thisStim) = CSC13_Wave(startplot:endplot);
    CSD_DataBlock (14,:,thisStim) = CSC14_Wave(startplot:endplot);
    CSD_DataBlock (15,:,thisStim) = CSC15_Wave(startplot:endplot);
    CSD_DataBlock (16,:,thisStim) = CSC16_Wave(startplot:endplot);
    CSD_DataBlock (17,:,thisStim) = CSC17_Wave(startplot:endplot);
    CSD_DataBlock (18,:,thisStim) = CSC18_Wave(startplot:endplot);
    CSD_DataBlock (19,:,thisStim) = CSC19_Wave(startplot:endplot);
    CSD_DataBlock (20,:,thisStim) = CSC20_Wave(startplot:endplot);
    CSD_DataBlock (21,:,thisStim) = CSC21_Wave(startplot:endplot);
    CSD_DataBlock (22,:,thisStim) = CSC22_Wave(startplot:endplot);
    CSD_DataBlock (23,:,thisStim) = CSC23_Wave(startplot:endplot);
    CSD_DataBlock (24,:,thisStim) = CSC24_Wave(startplot:endplot);
    CSD_DataBlock (25,:,thisStim) = CSC25_Wave(startplot:endplot);
    CSD_DataBlock (26,:,thisStim) = CSC26_Wave(startplot:endplot);
    CSD_DataBlock (27,:,thisStim) = CSC27_Wave(startplot:endplot);
    CSD_DataBlock (28,:,thisStim) = CSC28_Wave(startplot:endplot);
    CSD_DataBlock (29,:,thisStim) = CSC29_Wave(startplot:endplot);
    CSD_DataBlock (30,:,thisStim) = CSC30_Wave(startplot:endplot);
    CSD_DataBlock (31,:,thisStim) = CSC31_Wave(startplot:endplot);
    CSD_DataBlock (32,:,thisStim) = CSC32_Wave(startplot:endplot);

end

% Prep for offsetgraph for plot
offsetgraph = 500;

% Save CSD_Data
CSD_Data = mean(CSD_DataBlock,3);
save CSD_load.mat CSD_Data;