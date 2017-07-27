% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% Last edited: July 5, 2016

% NLX timestamp data is recorded in a 1D matrix
% For whatever reason, each Timestamp has 512 CSC_Samples associated with it

% This function reshapes this Timestamp data into a 1D matrix, which should
% correlate every CSC_Sample point with a discrete Timestamp


function timeVector = NB_timeVector3 (CSC_Timestamp)

% Creates 2D (512 x # of CSC_Timestamps) matrix with zeros
allTimeStamps = zeros(512,length(CSC_Timestamp));

% Loads CSC_Timestamp data into row 1 of the allTimeStamps matrix
allTimeStamps (1,:) = CSC_Timestamp;

% Calculates time in between each timestamp sampled
timeBetweenSamples = unique(diff (CSC_Timestamp));

% Number of points (CSC_Samples taken) in between each Timestamp
pointsPerSample = 512;

% Caclulates time in between each 512 points within a sample
timeBetweenPoints = timeBetweenSamples(1)./pointsPerSample;

for row = 2:512
    
    % Adds in Timestamp data for each CSC_Sample to the 2D Matrix
    % Adds at previous row's TimeStamp and timeBetweenPoints to current index
    allTimeStamps (row,:) = allTimeStamps(row-1,:) + timeBetweenPoints;
    
end

% Reshapes 2D matrix into 1D matrix; each timestamp in timeVector should be
% correlated with each point in CSC_Wave
timeVector = reshape (allTimeStamps, 1, 512*(length(CSC_Timestamp)));