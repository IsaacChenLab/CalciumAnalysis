% Nadir Bilici
% nadir.bilici@uphs.upenn.edu
% Last edited: July 7, 2016

% Plots the averaged responses from all 32 channels in msec time scale.
% Modify offsetgraph to get a better separation of the responses.

offsetgraph = 1000;
figure;
hold on

 for i =1:32;  
    
    % Normal CSD Plot of 32 channels
    x = ((1:length(CSD_Data(i,:)))/32)+offsetms;
    y = CSD_Data (i,:) - i*offsetgraph;
    plot(x,y);

    % Calculate slope of lines to determine PP facilitation/depression
    
    % Take difference between each two points of CSD_Data = first
    % derivative = slope between each two timepoints 
    a = diff(CSD_Data(i,:));
    % Smoothes the first derivative curve
    b = smooth(a,50);
    
%     figure;
%     plot(a);
%     plot(b);

    % Calculates minimum slope by determining the smallest value of b between b(startTime:endTime)
    slope1 = min(b(5:40));
    slope2 = min(b(55:90));
    
      ratio_down(i)=slope2/slope1;

    % Calculates maximum slope by determining the largest value of b between b(startTime:endTime)
    slope3 = max(b(5:40));
    slope4 = max(b(55:90));
    
      ratio_up(i) = slope4/slope3;
    
 end
 axis tight;

 down = transpose(ratio_down);
 up   = transpose(ratio_up);
 
 meanDown = mean(down);
 meanUp = mean(up);
 
save Ratio_Up.mat up;
save Ratio_Down.mat down;
 