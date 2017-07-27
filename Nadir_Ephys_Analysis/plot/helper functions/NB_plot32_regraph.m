% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% July 16, 2016

% Edited from plot32.m

% Plots the averaged responses from all 32 channels in msec time scale.
% Modify offsetgraph to get a better separation of the responses.

clf;

hold on

for i = 1:32;  
   
    % Calculates x coordinate in miliseconds
    x = ((1:length(CSD_Data(i,:)))/32)+offsetms;
    
    % Calculates y coordinate based on CSD_Data (which is based on averaged 
    % CSC_wave data
    y = CSD_Data (i,:) - i*offsetgraph;
    
    h = plot(x,y);
end

axis tight;