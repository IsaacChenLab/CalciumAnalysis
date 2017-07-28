function binMatrix = FR_vs_time(startTime, endTime, binSize, fileName)

%   Must be run in the same directory as a Spikes.mat file

% INPUT 
%   startTime,endTime = boundaires of time period to be analyzed (in seconds)
%   binSize = length of time for each bin
%   fileName = name of .mat file to which output will be saved (don't
%       include '.mat' at the end)

% OUTPUT
%   binMatrix = C x N matrix where c is the number of cells and N is the
%       number of complete bins in the time period ananlyzed. output(i,j)
%       is the aaverahe spikes/second of cell i throughout bin j. binMatrix
%       is both saved to [fileName].mat and returned by the function

load('Spikes.mat', 'dvSpikes')

startTime = startTime*1000000;
endTime = endTime*1000000;
binSize = binSize*1000000;

cellCount = size(dvSpikes.units, 2);
numBins = (endTime - startTime) / binSize;
binMatrix = zeros(cellCount, numBins);

for currBin = 1:numBins

    binStart = startTime + binSize*(currBin-1);
    binEnd = binStart + binSize;
    
    for c = 1:cellCount
        cSpikes = dvSpikes.units(c).stamps;
        spikesInBin = length( find( cSpikes>=binStart & cSpikes<binEnd));
        binMatrix(c,currBin) = (spikesInBin / binSize)*1000000;
    end
end

save( strcat(fileName,'.mat'), 'binMatrix');

end
