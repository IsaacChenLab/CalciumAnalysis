function PolarPlots(outputFolder, cellsToPlot, screenOff, scale, dataMatrix)

% INPUT
%   outputFolder = name (in SINGLE quotes) of output folder which will be 
%       created, and into which all of the output will be saved. If outputFolder 
%        = 'dont save' then no jpgs or .mat files will be saved and the script
%       will run much faster. If outputFolder is a complete path, then  '@'
%       should be added to the beginning so that user won't be prompted to
%       choose output directory later.
%   scale = the limit of the r axis for the plots (if the data exceeds r
%       for a given plot then the axis for that plot only will be extended)
%   cellsToPlot = a vector containing the cell number of each cell to be
%        plotted
%   screenOff = set to 'screenOff' if the columns for 'screenOff' are
%       included in the dataMatrix (ie the dataMatrix is N x 11 instead of N x
%       9, where N is the number of cells). The screenOff columns need to
%       be the first and last columns (these will be trimmed away). If
%       those columns are not present set screenOff to ''.
%   dataMatrix =  optional; n x 9 (or n x 11) matrix where n is # of cells.
%       The 9 columns are as follows: -180, -135, -90, -45, 0, 45, 90, 135,
%       180. Each entry is the firing rate in spikes/s. If dataMatrix is
%       left out then user will be prompted for a csv file containing the
%       matrix (csv should have no labels, just the matrix).

% OUTPUT
%   one polar graph saved as .jpg for each cell, saved to the folder that 
%   user provides


%get and read the csv
if ~exist('dataMatrix', 'var')  
    fprintf("\tSelect csv file...");
    [data_file, data_path] = uigetfile('*.csv','Select csv file...');
    fprintf("Selected!\n");
    M = csvread(strcat(data_path, data_file)); 
else
    M = dataMatrix;
end

if strcmp(screenOff, 'screenOff')
    M = M(2:end-1,:);
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


rads = deg2rad([-180,-135,-90,-45,0,45,90,135,180]);

for c = cellsToPlot
    
    %set up the figure
    name = sprintf('Cell_%.0f', c);
    f = figure('Name', name, 'NumberTitle','off');
    title(['Firing Rate vs Grating Orientation for Cell ' num2str(c)]);
    
    %set the r scale by creating white circle at that r
    polarplot(rads, scale*ones(size(rads)), '--w');  

    %plot and format
    hold on
    polarplot(rads, M(c,:));
    thetaticks(0:45:360);
    rticks([0 5 10 15 20 25]);
    hold off
    
    %save the file
    if ~strcmpi(outputFolder, 'dont save')
        saveas(f,strcat(target_folder, '/', name, '.jpg'));
    end
end

end
