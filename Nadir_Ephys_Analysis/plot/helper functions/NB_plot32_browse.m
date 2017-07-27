% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% Last edited: July 16, 2016

% Use this function to page through each stim

offsetgraph = 2000;

for thisStim = 1:numStims;

    for i = 1:32;  
        hold on
        
        % Calculates x coordinate in miliseconds
        x = ((1:lengthWave)/32)+offsetms;

        % Calculates y coordinate based on CSD_Data (which is based on averaged 
        % CSC_wave data
        y = CSD_DataBlock(i,:,thisStim) - i*offsetgraph;

        plot(x,y);
        
    end

    axis tight;
   
    hold off;
    

    % Ask user input to approve or reject plot    
    prompt = ['Currently looking at stim ', num2str(thisStim), '/', num2str(numStims), '. Press enter to go to next stim.']; 
    userInput = input(prompt);
    cla
    
end