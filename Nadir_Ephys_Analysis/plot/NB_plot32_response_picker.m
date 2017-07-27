% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% Last edited: July 11, 2016

% Edited from NB_plot32.m

% Plots the averaged responses from all 32 channels in msec time scale.
% Modify offsetgraph to get a better separation of the responses.

% After adding data to the block, plot this one temporarily

figure;
hold on;
offsetgraph = 1500;

for thisStim = 1:numStims;

    for i = 1:32;  
        hold on
        
        % Calculates x coordinate in miliseconds
        x = ((1:lengthWave)/32)+offsetms;

        % Calculates y coordinate based on CSD_Data (which is based on averaged 
        % CSC_wave data
        y = CSD_DataBlock(i,:,thisStim) - i*offsetgraph;

        tempPlot = plot(x,y);
        
    end

    axis tight;
   
    hold off;
    

    % Ask user input to approve or reject plot    
    prompt = ['Stim: ', num2str(thisStim), '/', num2str(numStims), '. Input 1 for positive, 0 for negative: ']; 
    userInput = input(prompt);
 
    % Copy elements of CSD_DataBlock to CSD_DataBlock_Positive or CSD_DataBlock_Negative
    for chan = 1:32;
        % If user approves
        if userInput == 1
            CSD_DataBlock_Positive (chan,:,thisStim) = CSD_DataBlock (chan,:,thisStim);
        else
            CSD_DataBlock_Negative (chan,:,thisStim) = CSD_DataBlock (chan,:,thisStim);
        end 
    end

    responseRate = (sum(rateHistogramArray))/thisStim;
    
    if userInput == 1
        rateHistogramArray(thisStim) = 1;
        disp (strcat('Positive. Current response rate: ', num2str(responseRate)))
    else 
        disp (strcat('Negative. Current response rate: ', num2str(responseRate)));
    end
    
    % Clear figure
    clf

end


CSD_Data_Positive = mean(CSD_DataBlock_Positive,3);
CSD_Data_Negative = mean(CSD_DataBlock_Negative,3);

save response_rate responseRate;
save CSD_load_response.mat rateHistogramArray;
save CSD_load_positive.mat CSD_Data_Positive;
save CSD_load_negative.mat CSD_Data_Negative;