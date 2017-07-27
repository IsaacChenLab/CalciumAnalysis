% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% Last edited: July 7, 2016

% Prepares the matrix of averaged responses after Paired Pulse (PP) protocol
% for all 32 channels. 

%Total period of time to be plotted
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

% Finds all timestamps when stim1 is turned on based on TTL events
stimOnTimestamps = EV_Timestamps(TTLon);

firstSetStimOnTimestamps = stimOnTimestamps(1:10);
secondSetStimOnTimestamps = stimOnTimestamps(11:20);

% Initialize vector to hold matched stim trials from first and second stim
% trials
matchedStimOnTimestamps = zeros(1,length(stimOnTimestamps));

for i = 1:2:19;
    
    % Matches stims 1,11; 2,12; 3,13 ... together into
    % matchedStimOnTimestamps. So, stims 1,11 are now index 1,2 in this vector
    setIndex = (round(i/2,0));
    
    matchedStimOnTimestamps(i) = firstSetStimOnTimestamps(setIndex);
    matchedStimOnTimestamps(i+1) = secondSetStimOnTimestamps(setIndex);    
    
end

% Number of stims calculated by number of times first stim was turned on
numstims = length(stimOnTimestamps);

% Interpolates the index of timeVector to match the timestamp of firstStimOnTimestamps
% interp_index is used to match firstStimOnTimestamps with CSC_Wave data
% This allows us to correlate stim events with CSC data
interp_index = interp1(timeVector,1:length(timeVector),matchedStimOnTimestamps,'nearest');

% Creates CSD_Datablock with zeros based on # channels and length of plot
CSD_DataBlock = zeros (channels, lengthWave, numstims);

for i = 1:length(matchedStimOnTimestamps);
    
    startstim = interp_index (i)+ offset;
    endstim = startstim + lengthWave - 1;

CSD_DataBlock (1,:,i) = CSC1_Wave(startstim:endstim);
CSD_DataBlock (2,:,i) = CSC2_Wave(startstim:endstim);
CSD_DataBlock (3,:,i) = CSC3_Wave(startstim:endstim);
CSD_DataBlock (4,:,i) = CSC4_Wave(startstim:endstim);
CSD_DataBlock (5,:,i) = CSC5_Wave(startstim:endstim);
CSD_DataBlock (6,:,i) = CSC6_Wave(startstim:endstim);
CSD_DataBlock (7,:,i) = CSC7_Wave(startstim:endstim);
CSD_DataBlock (8,:,i) = CSC8_Wave(startstim:endstim);
CSD_DataBlock (9,:,i) = CSC9_Wave(startstim:endstim);
CSD_DataBlock (10,:,i) = CSC10_Wave(startstim:endstim);
CSD_DataBlock (11,:,i) = CSC11_Wave(startstim:endstim);
CSD_DataBlock (12,:,i) = CSC12_Wave(startstim:endstim);
CSD_DataBlock (13,:,i) = CSC13_Wave(startstim:endstim);
CSD_DataBlock (14,:,i) = CSC14_Wave(startstim:endstim);
CSD_DataBlock (15,:,i) = CSC15_Wave(startstim:endstim);
CSD_DataBlock (16,:,i) = CSC16_Wave(startstim:endstim);
CSD_DataBlock (17,:,i) = CSC17_Wave(startstim:endstim);
CSD_DataBlock (18,:,i) = CSC18_Wave(startstim:endstim);
CSD_DataBlock (19,:,i) = CSC19_Wave(startstim:endstim);
CSD_DataBlock (20,:,i) = CSC20_Wave(startstim:endstim);
CSD_DataBlock (21,:,i) = CSC21_Wave(startstim:endstim);
CSD_DataBlock (22,:,i) = CSC22_Wave(startstim:endstim);
CSD_DataBlock (23,:,i) = CSC23_Wave(startstim:endstim);
CSD_DataBlock (24,:,i) = CSC24_Wave(startstim:endstim);
CSD_DataBlock (25,:,i) = CSC25_Wave(startstim:endstim);
CSD_DataBlock (26,:,i) = CSC26_Wave(startstim:endstim);
CSD_DataBlock (27,:,i) = CSC27_Wave(startstim:endstim);
CSD_DataBlock (28,:,i) = CSC28_Wave(startstim:endstim);
CSD_DataBlock (29,:,i) = CSC29_Wave(startstim:endstim);
CSD_DataBlock (30,:,i) = CSC30_Wave(startstim:endstim);
CSD_DataBlock (31,:,i) = CSC31_Wave(startstim:endstim);
CSD_DataBlock (32,:,i) = CSC32_Wave(startstim:endstim);

end

dt = (1/sampleRate)*1000;

% Set up empty CSD_DataBlock_IO
CSD_DataBlock_IO = zeros(channels, lengthWave, length(firstSetStimOnTimestamps));

% For every loop, add the averaged stim data from the first and second sets into CSD_DataBlock_IO
% After this loop, CSD_DataBlock_IO contains the averaged data for each IO stim step
for i = 1:length(firstSetStimOnTimestamps);
    adder = 0;
    
    for chan = 1:channels;
        IO1 = CSD_DataBlock(chan,:,i+adder);
        IO2 = CSD_DataBlock(chan,:,i+adder+1);
        averagedIO = (CSD_DataBlock(chan,:,i+adder) + CSD_DataBlock(chan,:,i+adder+1))/2;        
        CSD_DataBlock_IO(chan,:,i) = averagedIO;
    end
    adder = adder+1;
end

CSD_Data = mean (CSD_DataBlock,3);

save CSD_load.mat CSD_Data;
save CSD_IO_load.mat CSD_DataBlock_IO;
