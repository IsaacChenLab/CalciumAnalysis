function [AC_analysis] = ChenAutoCorr(maxLag, binSize, cellsToPlot, localMaxWidth, binMatrix)

% FUNCTION ARGUMENTS 
%   maxLag = max time(s) to offset the autocorr
%   binSize = length of time for each bin in file to be analyzed
%   cellsToPlot = a vector containing the cell number of each cell to be
%        plotted
%   localMaxWidth = an autocorrelation local maximum 'm' is defined as a value
%       which is larger than all the values 'localMaxWidth' seconds before 'm'
%       and after 'm'. The units for localMaxWidth are seconds. Usually
%       best around 0.3.
%   binMatrix = optional; if you want to pass in the data as an argument
%       (like if this function is being called within another script) then you won't be
%       prompted to select a data file.

% IN-FUNCTION PROMPTS
%   1. .mat file which is output from FC_vs_time() containg a single variable
%       'binMatrix' which is a 'C x binSize' matrix where C is the number
%       of cells in the associated recording
%   2. folder where all the figures generated will be saved

% PLOTS GENERATED
%   1. plot of bin'd firing rate over time
%   2. the auto-correlelogram
%   --- figures are displayed and saved in the output folder ---

% OUTPUT
%   AC_analysis = an array of struct, one struct for each cell. Each struct
%       has two fields: 'Time_Corr' and 'TimeOfMax_Period'.
%   Time_Corr simply has data that was plotted in autocorelograms. The first
%       column is x values (ie time offsets in seconds) and the second column is
%       the correlation coefficient for each time offset.
%   TimeOfMax_Period: Column 1 has the time offset (positive only) where
%       each local max was achieved (see 'localMaxWidth' for definition of local
%       max). Column 2 has the amount time between the corresponding local max
%       and the next one.


%if a binMatrix wasn't given as an argument, prompt the user for a file
if ~exist('binMatrix', 'var')
    fprintf('Select file to be analyzed...\n');
    [data_file, data_path] = uigetfile('*.mat', 'Select .mat file');
    fprintf('Selected!\n');
    load(strcat(data_path, data_file));
end

%prompt for file for output to be saved
fprintf('\nSelect folder where output files should be placed...\n');
target_folder = uigetdir('', 'Select output folder');
fprintf('Selected!\n');

%set some variables
numBins = size(binMatrix,2);
numLags = maxLag/binSize;
width = floor(localMaxWidth/binSize);
AC_analysis = cell(length(cellsToPlot),1);

for c = cellsToPlot
    
    %compute autocorrelelogram
    r_vector = xcorr(binMatrix(c,:), numLags, 'coeff');
    auto_x = -maxLag:binSize:maxLag;
    
    %plot firing rate
    fName = strcat('Activity_Cell_',num2str(c));
    f = figure('Name', fName, 'NumberTitle', 'off');
    ax1 = axes;
    plot(ax1, ((1:numBins)*binSize), binMatrix(c,:))
    xlabel(ax1, 'Time (s)');
    ylabel(ax1, 'Firing Rate (spikes/s)');
    
    %plot autocorelelogram
    aName = strcat('AutoCorr_Cell_',num2str(c));
    a = figure('Name', aName, 'NumberTitle', 'off');
    ax2 = axes;
    plot(ax2, auto_x, r_vector);
    xlabel(ax2, 'Time offset (s)');
    ylabel(ax2, 'Correlation Coeff');
    
    %save each figure
    saveas(f, strcat(target_folder, '/', fName, '.jpg'));
    saveas(a, strcat(target_folder, '/', aName, '.jpg'));
    
    %find the local maxima
    maxima = [];
    for x = (1+numLags):(length(r_vector)-width)
        if sum(r_vector(x) >= r_vector(x-width:x+width)) == width*2+1
            maxima = [maxima ; x];
        end
    end
    
    %find the periods between the local maxima
    numMax = length(maxima);
    periods = zeros(numMax,1);
    for i = 1:numMax-1
        periods(i) = maxima(i+1)-maxima(i);
    end
    
    %convert everything into column vectors, and convert units
    maxima = maxima - numLags - 1;
    maxima = maxima * binSize;
    periods = periods * binSize;
    r_vector = r_vector';
    auto_x = auto_x';
    
    %add the struct for this neuron to the array of structs
    s = struct('Time_Corr', [auto_x r_vector], 'TimeOfMax_Period', [maxima periods]);
    AC_analysis{c} = s;
end

end
