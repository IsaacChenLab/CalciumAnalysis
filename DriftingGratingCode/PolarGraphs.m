function PolarGraphs(scale)

% File contains n x 9 matrix where n is the number of cells
% The 9 columns are as follows: -180, -135, -90, -45, 0, 45, 90, 135, 180
% No labels in the csv
% Lone argument is max of r axis -- scale is same for every cell in the file

%get csv file
fprintf("\tSelect csv file...");
[data_file, data_path] = uigetfile('*.csv','Select csv file...');
fprintf("Selected!\n");

%get save folder
fprintf("\tSelect folder to save figures to...");
save_path = uigetdir('','Select folder to save figures to');
fprintf("Selected!\n");

M = csvread(strcat(data_path, data_file));
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