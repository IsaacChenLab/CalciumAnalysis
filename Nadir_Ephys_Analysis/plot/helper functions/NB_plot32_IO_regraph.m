% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% July 16, 2016

clf

hold on

for i = 1:10;  
   
    % Calculates x coordinate in miliseconds
    x = ((1:length(CSD_DataBlock_IO(32,:,i)))/32)+offsetms;
    
    % Calculates y coordinate based on CSD_Data (which is based on averaged 
    % CSC_wave data
    y = CSD_DataBlock_IO(32,:,i) - i*offsetgraph;
  
    plot1 = plot(x,y);

end

axis tight;