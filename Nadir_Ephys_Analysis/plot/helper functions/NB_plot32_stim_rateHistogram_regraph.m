% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% Last edited: July 18, 2016

% Clear positive and negative subplots
cla(subplot(10,10,[1 2 3 4 5 11 12 13 14 15 21 22 23 24 25 31 32 33 34 35 41 42 43 44 45 51 52 53 54 55 61 62 63 64 65 71 72 73 74 75 81 82 83 84 85]),'reset');
cla(subplot(10,10,[6 7 8 9 10 16 17 18 19 20 26 27 28 29 30 36 37 38 39 40 46 47 48 49 50 56 57 58 59 60 66 67 68 69 70 76 77 78 79 80 86 87 88 89 90]),'reset');

% Plot positive responses
for chan = 1:32;  
   
    % Calculates x coordinate in miliseconds
    x = ((1:length(CSD_Data(chan,:)))/32)+offsetms;
    
    % Calculates y coordinate based on CSD_Data (which is based on averaged 
    % CSC_wave data
    y = CSD_Data_Positive (chan,:) - chan*offsetgraph;
    
    hold on
    subplot(10,10,[1 2 3 4 5 11 12 13 14 15 21 22 23 24 25 31 32 33 34 35 41 42 43 44 45 51 52 53 54 55 61 62 63 64 65 71 72 73 74 75 81 82 83 84 85])    
    set(gca,'ytick',[])

    plot(x,y);
    title(['Positive (', num2str(sum(rateHistogramArray)),'/', num2str(length(rateHistogramArray)),' = ', num2str(responseRate) , ')' ]);
    axis tight
    hold off
end

for chan = 1:32;  
   
    % Calculates x coordinate in miliseconds
    x = ((1:length(CSD_Data(chan,:)))/32)+offsetms;
    
    % Calculates y coordinate based on CSD_Data (which is based on averaged 
    % CSC_wave data
    y = CSD_Data_Negative (chan,:) - chan*offsetgraph;
    
    hold on
    subplot(10,10,[6 7 8 9 10 16 17 18 19 20 26 27 28 29 30 36 37 38 39 40 46 47 48 49 50 56 57 58 59 60 66 67 68 69 70 76 77 78 79 80 86 87 88 89 90]);
    set(gca,'ytick',[])

    plot(x,y);
    title(['Negative (', num2str(length(rateHistogramArray) - sum(rateHistogramArray)),'/', num2str(length(rateHistogramArray)),' = ', num2str(1-responseRate) , ')' ]);
    axis tight
    hold off
end

