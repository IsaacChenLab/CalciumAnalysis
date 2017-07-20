function FiringRatePerEpoch(startTimes, fileName)

%   Must be run in the same directory as a Spikes.mat file

load('Spikes.mat', 'dvSpikes')

numEpochs = length(startTimes);
cellCount = size(dvSpikes.units, 2);

startTimesMicro = startTimes*1000000;
endTime = dvSpikes.stamps(end,1);

output = zeros(cellCount, numEpochs);

for c = 1:cellCount
    for x = 1:numEpochs
        
        spikes = dvSpikes.units(c).stamps;
        if x == numEpochs
            spikeCount = length( find( spikes > startTimesMicro(x)));
            FR = spikeCount / (endTime - startTimesMicro(x))*1000000;
        else
            spikeCount = length( find( spikes > startTimesMicro(x) & spikes<startTimesMicro(x+1)));
            FR = spikeCount / (startTimesMicro(x+1) - startTimesMicro(x))*1000000;
        end
        
        output(c,x) = FR;
    end
end

csvwrite(strcat(fileName,'.csv'), output);
%sxlswrite(strcat(fileName,'.xls'), output);

end