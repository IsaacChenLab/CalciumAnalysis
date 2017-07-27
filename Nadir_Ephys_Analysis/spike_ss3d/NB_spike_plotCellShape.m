% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% August 2016

% for each cluster, plot its spikes

numClusters = length(unique(TT1_CellNumbers));
figure % ('units','normalized','outerposition',[0 0 1 1]);
hold on

if numClusters > 1;
    % use subplots
    for cluster = 0:numClusters-1;
        for spike = 1:length(SamplesCSC1(1,:));
            if TT1_CellNumbers(spike) == cluster;
                subplot(1,numClusters,cluster+1);
                hold on
                plot(SamplesCSC1(1:32,spike));
            end
        end
        numSpikes = sum(TT1_CellNumbers(1,:) == cluster);
        title(['Cluster ' num2str(cluster) ' (n = ' num2str(numSpikes) ')']);   
    end
else
    % use single plot
    for spike = 1:length(SamplesCSC1(1,:));
        if TT1_CellNumbers(spike) == 0;
            hold on
            plot(SamplesCSC1(1:32,spike));
        end
        numSpikes = sum(TT1_CellNumbers(1,:) == 0);
        title(['Cluster ' num2str(0) ' (n = ' num2str(numSpikes) ')']);   
    end
end
