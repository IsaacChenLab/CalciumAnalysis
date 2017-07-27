% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% August 2016

% Adapted from MakeEvents.m
% Create timestamp arrays for different TTL signals.

% Find all TTL events when stim is turned on
TTLon = find (EV_TTLs==1);

% Finds all timestamps when stim is turned on based on TTL events
stimOnTimestamps = EV_Timestamps(TTLon);

% Number of stims calculated by number of times stim was turned on
numStims = length(stimOnTimestamps);