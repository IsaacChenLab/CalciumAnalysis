function ChenAutoCorr(cellMatrix, maxLag, binSize, cellsToPlot)

% INPUT 
%   cellMatrix
%        C x N matrix where c is the number of cells and N is the number of complete
%        bins in the time period ananlyzed. output(i,j) is the average spikes/second of
%        cell i throughout bin j
%   maxLag = max time(s) to offset the autocorr
%   binSize = length of time for each bin in cellMatrix
%   cellsToPlot = a vector containing the cell number of each cell to be
%        plotted

% OUTPUT for each cell
%   1. plot of firiring rate over time
%   2. the auto-correllogram

numBins = size(cellMatrix,2);
numLags = maxLag/binSize;

for i = cellsToPlot
    
    figure('Name', strcat("Activity of cell ",num2str(i)), 'NumberTitle', 'off');
    ax1 = axes;
    plot(ax1, ((1:numBins)*binSize), cellMatrix(i,:))
    xlabel(ax1, 'Time (s)');
    ylabel(ax1, 'Firing Rate (spikes/s)');
    
    figure('Name', strcat("AutoCorr of cell ",num2str(i)), 'NumberTitle', 'off');
    ax2 = axes;
    r = xcorr(cellMatrix(i,:), numLags, 'coeff');
    auto_x = -maxLag:binSize:maxLag;
    plot(ax2, auto_x, r);
    xlabel(ax2, 'Time offset (s)');
    ylabel(ax2, 'Correlation Coeff');
    
end

end