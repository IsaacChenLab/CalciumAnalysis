function output = FR_vs_time(startTime, endTime, binSize, cellCount)

% INPUT 
%   startTime, endTIme = boundaires of time period to be analyzed (in seconds)
%   binSize = length of time for each bin
%   cellCount = number of cells in file

% OUTPUT
%   C x N matrix where c is the number of cells and N is the number of complete
%   bins in the time period ananlyzed. output(i,j) is the aaverahe spikes/second of
%   cell i throughout bin j


load('Spikes.mat', 'dvSpikes')

startTime = startTime*1000000;
endTime = endTime*1000000;
binSize = binSize*1000000;

numBins = (endTime - startTime) / binSize;
output = zeros(cellCount, numBins);

for currBin = 1:numBins
   
    binStart = startTime + binSize*(currBin-1);
    binEnd = binStart + binSize;
    
    for c = 1:cellCount
        cSpikes = dvSpikes.units(c).stamps;
        spikesInBin = length( find( cSpikes>=binStart & cSpikes<binEnd));
        output(c,currBin) = (spikesInBin / binSize)*1000000;
    end
end

end
