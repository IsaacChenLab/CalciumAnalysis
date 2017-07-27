% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% Last edited: July 10, 2016

% Edited from plot32.m

% Modify offsetgraph to get a better separation of the responses.
offsetgraph = 1000;
targetChannel = 32;

plot1 = figure('units','normalized','outerposition',[0 0 1 1]);
hold on

for i = 1:10;  
   
    % Calculates x coordinate in miliseconds
    x = ((1:length(CSD_DataBlock_IO(targetChannel,:,i)))/32)+offsetms;
    
    % Calculates y coordinate based on CSD_Data (which is based on averaged 
    % CSC_wave data
    y = CSD_DataBlock_IO(targetChannel,:,i) - i*offsetgraph;
  
%    y = CSD_DataBlock_IO(32,:,i);
    plot1 = plot(x,y);

end

axis tight;
hold off;

%******************************************************************%
% Adjust these two values to adjust time for calculating amplitude %
msStartTime = 55;
msEndTime = 100;

plot2 = figure;
hold on

% For each step, find max value
for i = 1:10;
    ampX(i) = 50*i;
    yWave = CSD_DataBlock_IO(targetChannel,:,i);
    ampY(i) = max(yWave((msStartTime/1000*32000):(msEndTime/1000*32000)));
end

plot2 = plot(ampX,ampY,'b');

axis tight;

%axis([-50 500 (-10*10^4) 0]);
