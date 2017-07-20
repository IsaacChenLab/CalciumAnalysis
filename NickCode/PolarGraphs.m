function PolarGraphs(scale, dataMatrix)

% INPUT
%   scale = the maximum of the r axis for the graphs (if the data exceeds r
%       for a given graph then the axis for that graph only will be extended)
%   dataMatrix =  n x 9 matrix where n is the number of cells. The 9
%      columns are as follows: -180, -135, -90, -45, 0, 45, 90, 135, 180. Each
%      entry is the firing rate in spikes/s. dataMatrix does not need to be
%      given --- if left out then user will be prompted for a csv file
%      containing the matrix (no labels)

% OUTPUT
%   one polar graph saved as .jpg for each cell, saved to the file that the
%   user provides


if ~exist('dataMatrix', 'var')
    
    %get csv file
    fprintf("\tSelect csv file...");
    [data_file, data_path] = uigetfile('*.csv','Select csv file...');
    fprintf("Selected!\n");
    
    M = csvread(strcat(data_path, data_file));
    M = M(:,2:end-1);
else
    M = dataMatrix(:,2:end-1);
end

%get save folder
fprintf("\tSelect folder to save figures to...");
save_path = uigetdir('','Select folder to save figures to');
fprintf("Selected!\n");

cells = size(M,1);
rads = deg2rad([-180,-135,-90,-45,0,45,90,135,180]);

for c = 1:cells
    %set up the figure
    name = sprintf('Cell_%.0f', c);
    f = figure('Name', name, 'NumberTitle','off');
    title(name);
    
    %set the r scale by creating white circle at that r
    polarplot(rads, scale*ones(size(rads)), '--w');  

    %plot and format
    hold on
    polarplot(rads, M(c,:));
    thetaticks(0:45:360);
    rticks([0 5 10 15 20 25]);
    hold off
    
    %save the file
    saveas(f,strcat(save_path, '/', name, '.jpg'));
end

end
