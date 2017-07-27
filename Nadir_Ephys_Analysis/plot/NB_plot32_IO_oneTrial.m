% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% Last edited: July 10, 2016

% Edited from plot32.m

% Plots the averaged responses from all 32 channels in msec time scale.
% Modify offsetgraph to get a better separation of the responses.

offsetgraph = 3000;

plot1 = figure;

hold on

for i = 1:10;  
   
    % Calculates x coordinate in miliseconds
    x = ((1:length(CSD_DataBlock_IO(32,:,i)))/32)+offsetms;
    
    % Calculates y coordinate based on CSD_Data (which is based on averaged 
    % CSC_wave data
    y = CSD_DataBlock_IO(32,:,i) - i*offsetgraph;
  
%    y = CSD_DataBlock_IO(32,:,i);
    plot1 = plot(x,y);

end

axis tight;

NB_saveFigures;

hold off;

plot2 = figure;

hold on

% For each step, find max value
for i = 1:10;
    ampX(i) = 50*i;
    yWave = CSD_DataBlock_IO(32,:,i);
    ampY(i) = max(yWave((55/1000*32000):(100/1000*32000)));
end

plot2 = plot(ampX,ampY,'b');

axis tight;

%axis([-50 500 (-10*10^4) 0]);
