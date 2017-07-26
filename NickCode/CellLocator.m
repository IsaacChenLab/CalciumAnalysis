function Cell_Locations = CellLocator(outputFolder, cellCount, threshold, show_figs, selectivity)

% INPUT
%   outputFolder = name (in SINGLE quotes) of output folder which will be 
%       created, and into which all of the output will be saved. If outputFolder 
%        = 'dont save' then no .fig or .mat files will be saved and the script
%       will run much faster. If outputFolder is a complete path, then  '@'
%       should be added to the beginning so that user won't be prompted to
%       choose output directory later.
%   cellsCount = number of cells in the recording
%   threshold = optional; selectivity threshold that a cell must meet to be
%       considered "orientation selective". All cells have a selectivity metric
%       between 0 and 1 listed in Orientaton_Selectivity.mat. Default threshold
%       is 0.6
%   show_figs = optional; set to 'dont show' if you don't want the figure
%       windows to pop up. They'll still be saved unless 'dont save' is
%       specified
%  selectivity =  optional; the complete path for an
%       Orientation_Selectivity.mat file (output of PolarPlots.m). If not given,
%       user will be prompted to browse for the file.

% PLOTS GENERATED
%   1. Plot of relative locations of all the cells detected, based on which
%       probe they are nearest to. The orientation selective cells are red and
%       all other cells are blue.
%   --- figures are displayed and saved in the output folder ---

% OUTPUT
%   Cell_Locations is a struct with two fields: Selective_Cells and
%       NonSelective_Cells. Each field contains an matrix with 4 columns. Each
%       row represents a different cell.
%   Column 1 = the cell number
%   Column 2 = the channel (ie probe) to which that cell is closest
%   Column 3 = corresponding x coordinate
%   Column 4 = correspinding y coordinate

%handle some input
if exist('show_figs', 'var') && strcmp(show_figs, 'dont show')
   set(0,'DefaultFigureVisible','off');
else
    set(0,'DefaultFigureVisible','on');
end

if ~exist('threshold', 'var')
   threshold = 0.6;
end

% make sure that channelLocations.mat is in the same directory as this script
% this is where channels_xy comes from
load('channelLocations.mat');

%get and read the Orientation Selectivity file
if ~exist('selectivity', 'var')  
    fprintf("\tSelect orientation selectivity file...");
    [data_file, data_path] = uigetfile('*.mat','Select .mat file...');
    fprintf("Selected!\n");
    selectivity = strcat(data_path, data_file); 
end

%get and load the spikes.mat file
fprintf("\tSelect orientation selectivity file...");
[data_file, data_path] = uigetfile('*.mat','Select .mat file...');
fprintf("Selected!\n");
load(strcat(data_path, data_file));

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

load(selectivity);

%maxChannels = [9,11,6,24,9,30,24,18,6,27,13,26,27,12,17,5,6,31,11,4,5,18,20,22,20,27,11,25,13,17,15,10,9,28,8];

%set maxChannels and count repeats
maxChannels = zeros(cellCount, 1);
repeats = zeros(cellCount, 1);
for i = 1:cellCount
    maxChannels(i) = dvSpikes.units(i).maxChannels;
    for j = 1:i-1
        if maxChannels(j) == maxChannels(i)
            repeats(i) = repeats(i) + 1;
        end
    end
end

redCells = resultantVectors(resultantVectors(:,3) > threshold, 1);
allCells = 1:cellCount;
blueCells = setdiff(allCells, redCells);

red_repeats = repeats(redCells);
blue_repeats = repeats(blueCells);

redChannels = maxChannels(redCells);
blueChannels = maxChannels(blueCells);

red_xy = channels_xy(redChannels+1,2:3);
blue_xy = channels_xy(blueChannels+1,2:3);

% offset cells in the same location
red_repeats(red_xy(:,1) == 0) = red_repeats(red_xy(:,1) == 0) * -1;
blue_repeats(blue_xy(:,1) == 0) = blue_repeats(blue_xy(:,1) == 0) * -1;

red_xy_repeats = red_xy + [red_repeats*1.2 zeros(length(redCells),1)];
blue_xy_repeats = blue_xy + [blue_repeats*1.2 zeros(length(blueCells),1)];

%plot the points
f = figure('Name', 'Cell Locations', 'NumberTitle','off');
ax1 = axes;
hold(ax1, 'on');
scatter( ax1, red_xy_repeats(:,1), red_xy_repeats(:,2),25,... 
         'filled', 'r', 'DisplayName', 'Orientation Selective Cells');
scatter( ax1, blue_xy_repeats(:,1), blue_xy_repeats(:,2),25,...
         'filled', 'b',  'DisplayName', 'Non-Selective Cells');

% format
xlim(ax1, [-25, 50]);
ylim(ax1, [-200, 900]);
title(ax1, 'Relative Location of Detected Cells Based on Nearest Probe');
xlabel(ax1, 'Microns');
ylabel(ax1, 'Microns');
legend('show');

Cell_Locations = struct('Selective_Cells', [redCells redChannels' red_xy],...
                        'NonSelective_Cells', [blueCells' blueChannels' blue_xy]);

%save the output analysis
if ~strcmpi(outputFolder, 'dont save')
    saveas(f, strcat(target_folder, '/Cell_Locations.fig'));
    save(strcat(target_folder,'/Cell_Locations.mat'), 'Cell_Locations');
end

end

