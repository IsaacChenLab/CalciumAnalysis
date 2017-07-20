function ChenAutoCorr(maxLag, binSize, cellsToPlot)

% FUNCTION ARGUMENTS 
%   maxLag = max time(s) to offset the autocorr
%   binSize = length of time for each bin in file to be analyzed
%   cellsToPlot = a vector containing the cell number of each cell to be
%        plotted

% IN-FUNCTION PROMPTS
%   1. .mat file which is output from FC_vs_time() containg a single variable
%       'binMatrix' which is a 'C x binSize' matrix where C is the number
%       of cells in the associated recording
%   2. folder where all the figures generated will be saved

% OUTPUT for each cell
%   1. plot of bin'd firing rate over time
%   2. the auto-correlelogram
%   --- figures are displayed and saved in the output folder ---


fprintf('Select file to be analyzed...\n');
[data_file, data_path] = uigetfile('*.mat', 'Select .mat file');
fprintf('Selected!\n');

fprintf('\nSelect folder where output files should be placed...\n');
target_folder = uigetdir('', 'Select output folder');
fprintf('Selected!\n');

load(strcat(data_path, data_file));

numBins = size(binMatrix,2);
numLags = maxLag/binSize;

for i = cellsToPlot
    
    fName = strcat('Activity_Cell_',num2str(i));
    f = figure('Name', fName, 'NumberTitle', 'off');
    ax1 = axes;
    plot(ax1, ((1:numBins)*binSize), binMatrix(i,:))
    xlabel(ax1, 'Time (s)');
    ylabel(ax1, 'Firing Rate (spikes/s)');
    
    aName = strcat('AutoCorr_Cell_',num2str(i));
    a = figure('Name', aName, 'NumberTitle', 'off');
    ax2 = axes;
    r = xcorr(binMatrix(i,:), numLags, 'coeff');
    auto_x = -maxLag:binSize:maxLag;
    plot(ax2, auto_x, r);
    xlabel(ax2, 'Time offset (s)');
    ylabel(ax2, 'Correlation Coeff');
    
    saveas(f, strcat(target_folder, '/', fName, '.jpg'));
    saveas(a, strcat(target_folder, '/', aName, '.jpg'));
end

end
