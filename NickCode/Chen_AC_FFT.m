function AC_FFT_analysis = Chen_AC_FFT(outputFolder, startTime, endTime, binSize,...
                                       cellsToPlot, maxLag, localMaxWidth, binMatrixFile)

% FUNCTION ARGUMENTS
%   outputFolder = name (in SINGLE quotes) of output folder which will be 
%       created, and into which all of the output will be saved. If outputFolder 
%        = 'dont save' then no jpgs or .mat files will be saved and the script
%       will run much faster. If outputFolder is a complete path, then  '@'
%       should be added to the beginning so that user won't be prompted to
%       choose output directory later
%   startTime, endTime = define the period time (in seconds) over which to 
%       analyze. If endTime exceeds the final time point in the data file,
%       it is truncated appropriately
%   binSize = the amount of time associated with each bin in the file to be
%       analyzed
%   cellsToPlot = a vector containing the cell numbers of each cell to be
%        plotted
%   maxLag = max time(s) to offset the auto-correlation (recommended 5s)
%   localMaxWidth = an autocorrelation local maximum 'm' is defined as a
%       value which is larger than all the values 'localMaxWidth' seconds
%       before 'm' and after 'm'. The units for localMaxWidth are seconds
%       (recommended 3*binSize).
%   binMatrixFile = optional; if you want to pass in the complete path to
%       the data file as an argument (like if this function is being called
%       within another script) then you won't be prompted to select a data
%       file.

% IN-FUNCTION PROMPTS
%   1. .mat file which is output from FC_vs_time() containg a single variable
%       'binMatrix' which is a 'C x binSize' matrix where C is the number
%       of cells in the associated recording
%   2. folder where all the figures generated will be saved

% PLOTS GENERATED
%   1. Plot of bin'd firing rate over time
%   2. Auto-correlelogram
%   3. Fourier transform
%   --- figures are displayed and saved in the output folder ---

% OUTPUT
%   AC_FFT_Analysis = an array of struct, one struct for each cell. Each struct
%       has five fields:
%   AC_Time_Corr simply has data that was plotted in autocorelograms. Column 1
%       is x values (ie time offsets in seconds) and Column 2 is y values
%       (ie the correlation coefficient for each time offset).
%   AC_TimeOfMax_LocalMax_Period: Column 1 has the time offset (positive only) where
%       each local max was achieved (see 'localMaxWidth' for definition of local
%       max). Column 2 has the value of the correlation coeff which was
%       deemed a local max. Column 3 has the amount time between the
%       corresponding local max and the next one.
%   AC_ConfidenceInterval: upper and lower boundaries of the 95% confidence
%       intervals. If the distribution of spikes over time were truly random,
%       autocorrelation would be within the confidence interval (95% of the
%       time).
%   FFT_Freq_Amps simply has data that was plotted in FFT. Column 1
%       is x values (ie the frequencies) and Column 2 is y values
%       (ie the amplitude associated with each frequency).
%   FFT_MaxFreqs_MaxAmps: Column 1 has the 10 frequencies with greatest
%       amplitudes. Column 2 has the corresponding amplitudes.


%PROCESSING INPUT

%error if time input is invalid
if startTime < 0 || startTime > endTime
    error('invalid startTime given');
    return;
end

%if a binMatrix wasn't given as an argument, prompt the user for a file
if ~exist('binMatrixFile', 'var')
    fprintf('Select file to be analyzed...');
    [data_file, data_path] = uigetfile('*.mat', 'Select .mat file');
    fprintf('Selected!\n');
    load(strcat(data_path, data_file));
else
    load(binMatrixFile);
end

%prompt for file where output should be saved and create folder
if ~(strcmpi(outputFolder, 'dont save') || outputFolder(1) == '@')
    fprintf('Select folder where output files should be placed...');
    target_folder = uigetdir('', 'Select output folder');
    fprintf('Selected!\n');
    target_folder = strcat(target_folder, '/', outputFolder);
end

%if complete path was given as outputFolder argument, set target_folder to
%output folder
if outputFolder(1) == '@'
    target_folder = outputFolder(2:end);
end

%create output folder
if ~strcmpi(outputFolder, 'dont save')
   mkdir(target_folder);
end


%SETTING SOME VARIABLES

%shorten endTime if its too long
startBin = floor(startTime/binSize) + 1;
endBin = floor(endTime/binSize);
if endBin > size(binMatrix,2)
    endBin = size(binMatrix,2);
    endTime = endBin * binSize;
end

%set some general variables
time = endTime - startTime;
numBins = endBin - startBin;
numCells = length(cellsToPlot);
AC_FFT_analysis = cell(numCells,1);

%set some variables for autocorrelation
numLags = maxLag/binSize;
auto_x = -maxLag:binSize:maxLag;
width = floor(localMaxWidth/binSize);

%compute confidence intervals for autocorrelation
vcrit = sqrt(2)*erfinv(0.95);
lowCI = -vcrit/sqrt(numBins);
upCI = vcrit/sqrt(numBins);
lowCI_line = lowCI * ones(length(auto_x));
upCI_line = upCI * ones(length(auto_x));


for c = cellsToPlot
    
    data = binMatrix(c, startBin:endBin);
    
% FIRING RATE
    %plot firing rate
    FR_name = strcat('Activity_Cell_',num2str(c));
    FR_plot = figure('Name', FR_name, 'NumberTitle', 'off');
    ax1 = axes;
    plot(ax1, (startBin:endBin)*binSize, data)
    title(ax1, ['Firing Rate vs Time of Cell ', num2str(c)]);
    xlabel(ax1, 'Time (s)');
    ylabel(ax1, 'Firing Rate (spikes/s)');
    
    
% AUTOCORRELATION
    %compute autocorrelation
    r_vector = xcorr(data, numLags, 'coeff');
    
    %plot autocorrelation
    AC_name = strcat('AutoCorr_Cell_',num2str(c));
    AC_plot = figure('Name', AC_name, 'NumberTitle', 'off');
    ax2 = axes;
    plot(ax2, auto_x, r_vector);
    title(ax2, ['Autocorrelogram for Cell ', num2str(c)]);
    xlabel(ax2, 'Time offset (s)');
    ylabel(ax2, 'Correlation Coeff');
    
    %plot AC confidence intervals lines on the graph
    hold on
    plot(ax2, auto_x, lowCI_line,'r');
    plot(ax2, auto_x, upCI_line,'r');
    ylim([lowCI-0.03 1]);
    hold off
      
    %find AC local maxima
    maxima = [];
    maximaTimes = [];
    for x = (1+numLags):(length(r_vector)-width)
        if sum(r_vector(x) >= r_vector(x-width:x+width)) == width*2+1
            maxima = [maxima ; r_vector(x)]; 
            maximaTimes = [maximaTimes ; auto_x(x)];
        end
    end
    
    %find the periods between the local maxima
    numMax = length(maxima);
    periods = zeros(numMax,1);
    for i = 1:numMax-1
        periods(i) = maximaTimes(i+1) - maximaTimes(i);
    end
    
    
%FOURIER TRANSFORM
    %compute and normalize the single spectrum FFT
    F = fft(data);
    F = abs(F) / numBins;
    F = F(1: floor(numBins/2)+1);
    F(2:end) = F(2:end)*2;
    
    %set the X axis
    freqs = (0:numBins/2) / time;
    
    %plot the FFT
    FFT_name = strcat('FFT_Cell_',num2str(c));
    FFT_plot = figure('Name', FFT_name, 'NumberTitle', 'off');
    ax3 = axes;
    plot(ax3, freqs, F);
    title(ax3, ['Fourier Transform of Firing Rate vs Time for Cell ' num2str(c)]);
    xlabel(ax3, 'Frequency (Hz)');
    ylabel(ax3, 'Amplitude');
       
    %find the 10 strongest frequencies and corresponding amplitudes
    [ampRankValues, ampRankIndex] = sort(F,'descend');
    maxAmps = ampRankValues(1:10);
    maxFreqs = freqs(ampRankIndex(1:10));
    
    
% OUTPUT
    %save each figure
    if ~strcmpi(outputFolder, 'dont save')
        saveas(FR_plot, strcat(target_folder, '/', FR_name, '.jpg'));
        saveas(AC_plot, strcat(target_folder, '/', AC_name, '.jpg'));
        saveas(FFT_plot, strcat(target_folder, '/', FFT_name, '.jpg')); 
    end
    
    %add the struct for this neuron to the array of structs
    s = struct('AC_Time_Corr', [auto_x' r_vector'],...
               'AC_TimeOfMax_LocalMax_Period', [maximaTimes maxima periods],...
               'AC_ConfidenceIntervals', [lowCI; upCI],...
               'FFT_Freqs_Amps', [freqs' F'],...
               'FFT_MaxFreqs_MaxAmps', [maxFreqs' maxAmps']);
    AC_FFT_analysis{c} = s;        
end

%save the output analysis
if ~strcmpi(outputFolder, 'dont save')
    save( strcat(target_folder,'/','AC_FFT_analysis.mat'), 'AC_FFT_analysis');
end

end
