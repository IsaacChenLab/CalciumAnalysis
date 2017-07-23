function FFT_analysis = ChenFFT(outputFolder, startTime, endTime, binSize, cellsToPlot, binMatrix)


%error if time input is invalid
if startTime < 0 || startTime > endTime
    error('invalid startTime given');
    return;
end


%if a binMatrix wasn't given as an argument, prompt the user for a file
if ~exist('binMatrix', 'var')
    fprintf('Select file to be analyzed...');
    [data_file, data_path] = uigetfile('*.mat', 'Select .mat file');
    fprintf('Selected!\n');
    load(strcat(data_path, data_file));
end

%prompt for file for output to be saved
fprintf('\nSelect folder where output files should be placed...');
target_folder = uigetdir('', 'Select output folder');
fprintf('Selected!\n');

%create the output folder
target_folder = strcat(target_folder, '/', outputFolder);
mkdir(target_folder);


%shorten endTime if its too long
startBin = floor(startTime/binSize) + 1;
endBin = floor(endTime/binSize);
if endBin > size(binMatrix,2)
    endBin = size(binMatrix,2);
    endTime = endBin * binSize;
end

%set some variables
time = endTime - startTime;
numBins = endBin - startBin;
numCells = length(cellsToPlot);
FFT_analysis = cell(numCells,1);


for c = cellsToPlot

    %compute and normalize the single spectrum FFT
    data = binMatrix(c, startBin:endBin);
    F = fft(data);
    F = abs(F) / numBins;
    F = F(1: floor(numBins/2)+1);
    F(2:end) = F(2:end)*2;
    
    %set the X axis
    freqs = (0:numBins/2) / time;
    
    %plot the FFT
    name = strcat('FFT_Cell',num2str(c));
    f = figure('Name', name, 'NumberTitle', 'off');
    plot(freqs, F);
    title(['Fourier Transform of Firing Rate vs Time for Cell ' num2str(c)]);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    
    %save each figure
    saveas(f, strcat(target_folder, '/', name, '.jpg'));
    
    %find the 10 strongest frequencies and corresponding amplitudes
    [ampRankValues, ampRankIndex] = sort(F,'descend');
    maxAmps = ampRankValues(1:10);
    maxFreqs = freqs(ampRankIndex(1:10));
    
    %add the struct for this neuron to the array of structs
    s = struct('Freqs_Amps', [freqs' F'],...
               'MaxFreqs_MaxAmps', [maxFreqs' maxAmps']);
    FFT_analysis{c} = s;
    
end

    %save the array of structs
    save( strcat(target_folder,'/','FFT_analysis.mat'), 'FFT_analysis');

end
