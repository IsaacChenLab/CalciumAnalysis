function CellLocator(cellCount, thresh, selectivity)

%make sure that channelLocations.mat is in the same directory as this
load('channelLocations.mat');

%get and read the csv
if ~exist('selectivity', 'var')  
    fprintf("\tSelect orientation selectivity file...");
    [data_file, data_path] = uigetfile('*.mat','Select .mat file...');
    fprintf("Selected!\n");
    selectivity = strcat(data_path, data_file); 
end

load(selectivity);
maxChannels = [9,11,6,24,9,30,24,18,6,27,13,26,27,12,17,5,6,31,11,4,5,18,20,22,20,27,11,25,13,17,15,10,9,28,8];

repeats = zeros(cellCount, 1);
for i = 1:cellCount
    for j = 1:i-1
        if maxChannels(j) == maxChannels(i)
            repeats(i) = repeats(i) + 1;
        end
    end
end


redCells = resultantVectors(resultantVectors(:,3) > thresh, 1);
allCells = 1:cellCount;
blueCells = setdiff(allCells, redCells);

red_repeats = repeats(redCells);
blue_repeats = repeats(blueCells);

%redChannels = dvSpikes.units.maxChannels(redCells) + 1;
%blueChannels = dvSpikes.units.maxChannels(blueCells) + 1;

redChannels = maxChannels(redCells) + 1;
blueChannels = maxChannels(blueCells) + 1;

red_xy = channels_xy(redChannels,2:3);
blue_xy = channels_xy(blueChannels,2:3);

red_repeats(red_xy(:,1) == 0) = red_repeats(red_xy(:,1) == 0) * -1;
blue_repeats(blue_xy(:,1) == 0) = blue_repeats(blue_xy(:,1) == 0) * -1;

red_xy_repeats = red_xy + [red_repeats*1.2 zeros(length(redCells),1)];
blue_xy_repeats = blue_xy + [blue_repeats*1.2 zeros(length(blueCells),1)];

%plot the points
f = figure('Name', 'Cell Locations', 'NumberTitle','off');
ax1 = axes;
hold(ax1, 'on');
scatter(ax1, red_xy_repeats(:,1), red_xy_repeats(:,2),25, 'filled', 'r', 'DisplayName', 'Orientation Selective Cells');
scatter(ax1, blue_xy_repeats(:,1), blue_xy_repeats(:,2),25, 'filled', 'b',  'DisplayName', 'Non-Selective Cells');

% format
xlim(ax1, [-25, 50]);
ylim(ax1, [-200, 900]);
title(ax1, 'Relative Location of Detected Cells Based on Nearest Probe');
xlabel(ax1, 'Microns');
ylabel(ax1, 'Microns');
legend('show');

end

